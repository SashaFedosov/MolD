"""
This script compiles sDNC-based DNA diagnoses for a pre-defined taxa in a dataset
"""
import os, sys
import random
import argparse
######################################################################## FUNCTIONS
#***STEP 1 - SORTING ENTRIES BY CLADE AND IDENTIFYING NUCLEOTIDE POSITIONS SHARED WITHIN CLADE
def Step1(raw_records):
    Clades=[]
    for i in range(len(raw_records)):
        Clade=raw_records[i][1]
        if Clade not in Clades:
            Clades.append(Clade)
    clade_sorted_seqs = {}
    for letter in Clades:
        clade_sorted_seqs[letter]=[]
        for i in range(len(raw_records)):
            if raw_records[i][1]==letter:
                clade_sorted_seqs[letter].append(raw_records[i][2])
    shared_positions={}
    for key in clade_sorted_seqs:
        sh_pos=[]
        for i in range(len(clade_sorted_seqs[key][0])):
            shared_nucleotide = True
            csm = clade_sorted_seqs[key][0][i] #candidate shared nucleotide
            for j in range(1, len(clade_sorted_seqs[key])):
                if clade_sorted_seqs[key][j][i] != csm:
                    shared_nucleotide = False
                    break
            if shared_nucleotide == True and csm != 'N':
                sh_pos.append(i)
        shared_positions[key]=sh_pos
    return Clades, clade_sorted_seqs, shared_positions

#***STEP 2 COMPILING COMPARISON LISTS FOR CLADES AND IDENTIFYING VARIABLE POSITIONS AND N PRIORITY POSITIONS WITH LARGEST CUTOFFS
def C_VP_PP(clade_sorted_seqs, clade, shared_positions, CUTOFF):# complist_variable_positions_priority_positions; Arguments: dictionary, string, dictionary
    CShN={}#a dictionary keys - clade shared positions, values - nucleotides at those positions
    for pos in shared_positions[clade]:
        CShN[pos] = clade_sorted_seqs[clade][0][pos]#creates a dictionary shared position : nucleotide
    complist=[]
    for key in clade_sorted_seqs:
        if key != clade:
            complist = complist + clade_sorted_seqs[key]#creates a list of all other sequences for comparison
    cutoffs = {}
    pures = []####
    for key in CShN:
        newcomplist = []
        for k in complist:
            if k[key] == CShN[key]:
                newcomplist.append(k)
            else: continue
        cutoffs[key] = len(complist) - len(newcomplist)
        if len(newcomplist) == 0:####
            pures.append(key)####
    CPP = []
    for key in sorted(cutoffs, key = cutoffs.get, reverse = True):
        CPP.append(key)
    Clade_priority_positions = {}
    for position in CPP[:CUTOFF]:#Here you define how many of the clade shared combinations are used in subsequent search
        Clade_priority_positions[position] = CShN[position]
    return complist, Clade_priority_positions, cutoffs, pures####

#***STEPS 3 RANDOM SEARCH ACROSS PRIORITY POSITIONS TO FIND RAW DIAGNOSTIC COMBINATIONS AND TO SUBSEQUENTLY REFINE THEM
def random_position(somelist, checklist):#gives a random index (integer) of the specified range, and returns indexed somelist element if it is not present in the checklist 
    while True:
        i = random.randint(0, len(somelist) - 1)
        if somelist[i] not in checklist:
            return somelist[i]
            break
        else:
            continue

def step_reduction_complist(clade, complist, CPP, checked_ind):#checks randomly selected positions of CladeSharedNucleotides with sequences of other clades, until a diagnostic combination of nucleotides for a selected clade is found.
    if len(complist) == 0:
        return checked_ind
    elif len(checked_ind) == len(CPP):
        return checked_ind
    else:
        newcomplist = []
        pos = random_position(list(CPP.keys()), checked_ind)
        for j in complist:
            if j[pos] == CPP[pos]:
                newcomplist.append(j)
            else: continue
        new_checked_ind = checked_ind + [pos]
        return step_reduction_complist(clade, newcomplist, CPP, new_checked_ind)

def ConditionD(newcomb, complist, CPP):#The function checks the 'Condition D' - i.e. whither any given combination of nucleotide positions is diagnostic for the selected clade
    ContD = False
    for i in newcomb:
        newcomplist = []
        for m in complist:
            if m[i] == CPP[i]:
                newcomplist.append(m)
            else: continue
        complist = newcomplist
    if len(complist) == 0:
        ContD = True
    return ContD

def RemoveRedundantPositions(raw_comb, complist, CPP):# The function removes positions from the raw combinations one by one, and then checks whether new combination fulfills the condition D, thus recursively reducing the diagnostic combination.
    red_possible = False
    for j in raw_comb:
        newcomb = [k for k in raw_comb if k != j]
        if ConditionD(newcomb, complist, CPP) == True:
            red_possible = True
            return RemoveRedundantPositions(newcomb, complist, CPP)
        else: pass
    if red_possible == False:
        return raw_comb

#PUTS EVERYTHING TOGETHER - 20000 ROUNDS OF RANDOM SEARCH FOLLOWED BY REFINING OF 500 SHORTEST COMBINATIONS
def Diagnostic_combinations(qCLADE, complist, CPP, n1, maxlen1, maxlen2):
    Achecked_ind = []
    bestlists = []
    n = n1
    while n>0:#STEP3 proposes raw diagnostic combinations
        m = step_reduction_complist(qCLADE, complist, CPP, Achecked_ind)
        if len(m) < maxlen1 and sorted(m) not in bestlists:
            bestlists.append(sorted(m))
        n=n-1
    bestlists.sort(key=len)
    priority_lists = bestlists[:500]
    diagnostic_combinations = []
    for raw_comb in priority_lists:
        refined_comb = RemoveRedundantPositions(raw_comb, complist, CPP)
        if not refined_comb in diagnostic_combinations and len(refined_comb) <=maxlen2:
            diagnostic_combinations.append(refined_comb)
    diagnostic_combinations.sort(key=len)
    return diagnostic_combinations

#***STEP 4 ANALYSIS OF OUTPUT sDNCs
def IndependentKey(diagnostic_combinations):#PRESENTLY NOT INVOLVED - returns independent diagnostic nucleotide combinations, and identifies key nucleotide positions
    independent_combinations = []
    selected_positions = []
    for i in range(len(diagnostic_combinations)):
        if len(selected_positions) == 0:
            for j in range(0, i):
                if len(set(diagnostic_combinations[i]) & set(diagnostic_combinations[j])) == 0 and len(set(diagnostic_combinations[i]) & set(selected_positions)) == 0:
                    independent_combinations.append(diagnostic_combinations[i])
                    independent_combinations.append(diagnostic_combinations[j])
                    for k in range(len(diagnostic_combinations[i])):
                        selected_positions.append(diagnostic_combinations[i][k])
                    for l in range(len(diagnostic_combinations[j])):
                        selected_positions.append(diagnostic_combinations[j][l])
        else:
            if len(set(diagnostic_combinations[i]) & set(selected_positions)) == 0:
                independent_combinations.append(diagnostic_combinations[i])
                for k in range(len(diagnostic_combinations[i])):
                    selected_positions.append(diagnostic_combinations[i][k])
    independent_combinations.sort(key=len)
    key_positions = []
    for pos in diagnostic_combinations[0]:
        KP = True
        for combination in diagnostic_combinations[1:]:
            if pos not in combination:
                KP = False
                break
            else: continue
        if KP == True:
            key_positions.append(pos)
    return independent_combinations, key_positions

#SPECIFIC FUNCTIONS FOR THE sDNCs
def random_sequence2(SEQ, Pdiff, key_pos):#Returns a string (DNA sequence) that is Pdiff% different from the input string. Argument: string, integer (%), list. This and subsequent functions are only used for sDNC-based diagnoses
    N = len(SEQ)*Pdiff/100
    PosToChange = random.sample(range(0, len(SEQ)), N)
    Nucleotides = ['A', 'C', 'G', 'T']
    NEWSEQ = ''
    for i in range(len(SEQ)):
        if i not in PosToChange or i in key_pos:
            NEWSEQ = NEWSEQ + SEQ[i]
        else:
            NEWSEQ = NEWSEQ + random.sample([j for j in Nucleotides if j != SEQ[i]], 1)[0]
    return NEWSEQ

def GenerateBarcode2(Diagnostic_combinations, length):#This function calculates diagnostic combinations and assembles a barcode of desired length for a query taxon. First all single position DNCs are added, then based on the frequency of a nucleotide position in the DNCs of the 2 positions, and then based on the frequency of a position in longer DNCs
    OnePos = []
    positions_involved2 = []
    positions_involved_more = []
    for comb in Diagnostic_combinations:
        if len(comb) == 1:
            OnePos.append(comb[0])
        elif len(comb) == 2:
            for i in comb:
                positions_involved2.append(i)
        else:
            for j in comb:
                positions_involved_more.append(j)
    Setin = []
    for pos in sorted(positions_involved2, key=positions_involved2.count, reverse = True):
        if not pos in Setin:
            Setin.append(pos)
    for pos1 in sorted(positions_involved_more, key=positions_involved_more.count, reverse = True):
        if not pos1 in Setin:
            Setin.append(pos1)
    Ordered = OnePos + Setin
    return Ordered[:length]

def Screwed_dataset31(dataset, Taxon, sp_per_clade_to_screw, Prop_to_screw, key_pos, Percent_difference, Cutoff):#implements a random_sequence function to build a random sequence dataset where 10% sequences of the dataset (but not more than 20 species per clade) are Pdiff percent different from the original one. Arguments: list, string, list, list.
    Pdiff_records=[]
    screwed_clades = {}
    seq_to_screw = random.sample(range(len(dataset)), int(len(dataset)*Prop_to_screw))
    for i in range(len(dataset)):
        if i in seq_to_screw:
            clade_record = dataset[i][1]
            if clade_record not in list(screwed_clades.keys()):
                screwed_clades[clade_record] = 1                       
                Pdiff_record = [dataset[i][0], dataset[i][1], random_sequence2(dataset[i][2], Percent_difference, key_pos)]
                Pdiff_records.append(Pdiff_record)
            else:
                if screwed_clades[clade_record] < sp_per_clade_to_screw:
                    screwed_clades[clade_record] += 1                       
                    Pdiff_record = [dataset[i][0], dataset[i][1], random_sequence2(dataset[i][2], Percent_difference, key_pos)]
                    Pdiff_records.append(Pdiff_record)
                else:
                    Pdiff_records.append(dataset[i])
        else:
            Pdiff_records.append(dataset[i])
                
    Clades, clade_sorted_seqs, shared_positions = Step1(Pdiff_records)
    x,y,z,pures = C_VP_PP(clade_sorted_seqs, Taxon, shared_positions, Cutoff)#STEP2####
    return x, y

################################################READ IN PARAMETER FILE AND DATA FILE

def get_args(): #arguments needed to give to this script
    parser = argparse.ArgumentParser(description="run MolD")
    required = parser.add_argument_group("required arguments")
    required.add_argument("-i", help="textfile with parameters of the analysis", required=True)
    return parser.parse_args()

def main():
    args = get_args()
    ParDict = {}
    with open(args.i) as params:
        for line in params:
            line = line.strip()
            if line.startswith('#'):
                pass
            else:
                if len(line.split('=')) == 2 and len(line.split('=')[1]) != 0:
                    ParDict[line.split('=')[0]] = line.split('=')[1]

    f = open(ParDict['INPUT_FILE'], 'r') 
    imported=[]#set up a new dictionary with species and identifiers
    for line in f:
        line=line.rstrip()
        words=line.split()
        if len(words) != 3:
            print 'Check number of entries in', words[0]
            #break
        else:
            imported.append([words[0], words[1], words[2].upper()])
    f.close()
    if len(set([len(i[2]) for i in imported])) != 1:
        print 'Alignment contains sequences of different lengths:', set([len(i[2]) for i in imported])
    else:
        FragmentLen = len(imported[0][2])
    if 'NumberN' in list(ParDict.keys()):#How many ambiguously called nucleotides are allowed
        NumberN = int(ParDict['NumberN'])
    else:
        NumberN = 5
    raw_records=[]
    for i in imported:
        if i[2].count('N') < NumberN and len(i[2]) == FragmentLen:
            raw_records.append(i)
    print 'Maximum undetermined nucleotides allowed:', NumberN
    print 'Length of the alignment:', FragmentLen    
    print 'Read in', len(raw_records), 'sequences'
   
    #############################################READ IN OTHER ANALYSIS PARAMETERS
    if ParDict['qTAXA'] in ['ALL', 'All', 'all']:#qTAXA
        qCLADEs = []
        for i in raw_records:
            if not i[1] in qCLADEs:
                qCLADEs.append(i[1])
    elif ParDict['qTAXA'][0] == '>':
        NumSeq = int(ParDict['qTAXA'][1:])
        Taxarecords = [i[1] for i in raw_records]
        qCLADEs = []
        for j in Taxarecords:
            if Taxarecords.count(j) >= NumSeq and not j in qCLADEs:
                qCLADEs.append(j)
    else:
        qCLADEs = ParDict['qTAXA'].split(',')
    print 'focus taxa:', qCLADEs, len(qCLADEs)

    if 'Cutoff' in list(ParDict.keys()):#CUTOFF Number of the informative positions to be considered, default 100
        Cutoff = int(ParDict['Cutoff'])
    else:
        Cutoff = 100
    print 'Cutoff set as:', Cutoff
    if 'Number_of_iterations' in list(ParDict.keys()):#Number iterations of MolD
        N1 = int(ParDict['Number_of_iterations'])
    else:
        N1 = 10000
    print 'Number iterations of MolD set as:', N1
    
    if 'MaxLen1' in list(ParDict.keys()):#Maximum length for the raw pDNCs
        MaxLen1 = int(ParDict['MaxLen1'])
    else:
        MaxLen1 = 12
    print 'Maximum length of raw pDNCs set as:', MaxLen1
    
    if 'MaxLen2' in list(ParDict.keys()):#Maximum length for the refined pDNCs
        MaxLen2 = int(ParDict['MaxLen2'])
    else:
        MaxLen2 = 7
    print 'Maximum length of refined pDNCs set as:', MaxLen2
    
    if int(ParDict['Taxon_rank']) == 1:#read in taxon rank to configure Pdiff parameter of artificial dataset
        Percent_difference = 1
    else:
        Percent_difference = 3
    print 'artificial sequences are', Percent_difference, 'diverged from original ones'
    
    if 'PrSeq' in list(ParDict.keys()):#Proportion of sequences in the dataset to be modified for the artificial dataset construction
        P_to_screw = float(ParDict['PrSeq'])
    else:
        if Percent_difference == 1 and len(raw_records) < 500:
            P_to_screw = 0.5
        elif Percent_difference == 1 and len(raw_records) >= 500:
            P_to_screw = 0.25
        else:
            P_to_screw = 0.1
    print 'Proportion of modified sequences in artificial dataset', P_to_screw
    
    if 'NMaxSeq' in list(ParDict.keys()):#Maximum number of sequences per taxon to be modified
        Seq_per_clade_to_screw = int(ParDict['NMaxSeq'])
    else:
        Seq_per_clade_to_screw = 10
    print 'Maximum number of modified sequences per clade', Seq_per_clade_to_screw
    
    if 'Scoring' in list(ParDict.keys()):
        if ParDict['Scoring'] == 'lousy':
            threshold = 66
        elif ParDict['Scoring'] == 'moderate':
            threshold = 75
        elif ParDict['Scoring'] == 'stringent':
            threshold = 90
        elif ParDict['Scoring'] == 'very_stringent':
            threshold = 95
        else:
            threshold = 75
    else:
        threshold = 75
    print 'Scoring of the sDNCs:', ParDict['Scoring'], 'Threshold in two consequtive runs:', threshold
    
    ###################################################IMPLEMENTATION
    #Setting up a new class just for the convenient output formatting
    class SortedDisplayDict(dict):#this is only to get a likable formatting of the barcode
        def __str__(self):
            return "[" + ", ".join("%r: %r" % (key, self[key]) for key in sorted(self)) + "]"
    
    #Calling functions and outputing results
    g = open(ParDict['OUTPUT_FILE'], "w")#Initiating output file
    for qCLADE in qCLADEs:
        print '**************', qCLADE, '**************'
        print >>g, '**************', qCLADE, '**************'
        Clades, clade_sorted_seqs, shared_positions = Step1(raw_records)#STEP1
        x,y,z,pures = C_VP_PP(clade_sorted_seqs, qCLADE, shared_positions, Cutoff)#STEP2 ####
        newy = {key:y[key] for key in y if not key in pures} ####
        print 'Sequences analyzed:', len(clade_sorted_seqs[qCLADE])
        print >>g, 'Sequences analyzed:', len(clade_sorted_seqs[qCLADE])
        ND_combinations = [[item] for item in pures] ####
        print 'single nucleotide pDNCs:', str(ND_combinations).replace('[', '').replace(']', '')
        N = 1
        while N > 0:#STEP3
            try:
                q = Diagnostic_combinations(qCLADE, x, newy, N1, MaxLen1, MaxLen2) ####
            except IndexError:
                print N, 'IndexError'
                continue
            for comb in q:
                if not comb in ND_combinations:
                    ND_combinations.append(comb)
            N-=1
        ND_combinations.sort(key=len)
        #################################### pDNC output
        try:
            Nind, KeyPos = IndependentKey(ND_combinations)#STEP4
        except IndexError:
            print qCLADE, 'no pDNCs recovered for', qCLADE
            continue
        Allpos = []#Create list of all positions involved in pDNCs
        for comb in ND_combinations:
            for pos in comb:
                if not pos in Allpos:
                    Allpos.append(pos)
        print 'pDNCs retrieved :', len(ND_combinations), 'Positions involved:', len(Allpos), 'Independent pDNCs:', len(Nind)
        print >>g, 'pDNCs retrieved :', len(ND_combinations), 'Positions involved:', len(Allpos), 'Independent pDNCs:', len(Nind)
        print 'Shortest retrieved diagnostic combination:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in ND_combinations[0]]})
        print >>g, 'Shortest retrieved diagnostic combination:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in ND_combinations[0]]})
        ######################################################## sDNC output
        Barcode_scores = []#Initiate a list for sDNC scores
        npos = len(ND_combinations[0])
        while npos <= min([10, len(Allpos)]):#in this loop the positions are added one-by-one to a sDNC and the sDNC is then rated on the artificially generated datasets
            Barcode = GenerateBarcode2(ND_combinations, npos)#Initiate a sDNC
            Barcode_score = 0#Initiate a score to rate a sDNC
            N = 100
            while N > 0:
                NComplist, NCPP = Screwed_dataset31(raw_records, qCLADE, Seq_per_clade_to_screw, P_to_screw, KeyPos, Percent_difference, Cutoff)#Create an artificial dataset
                NBarcode = [i for i in Barcode if i in list(NCPP.keys())]
                if len(Barcode) - len(NBarcode) <= 1 and ConditionD(NBarcode, NComplist, NCPP) == True:#Check whether the sDNC works for the artificial dataset NEW
                    Barcode_score +=1
                N -=1
            print npos, 'sDNC_score (100):', Barcode, '-', Barcode_score
            print >>g, npos, 'sDNC_score (100):', Barcode, '-', Barcode_score
            Barcode_scores.append(Barcode_score)
            if len(Barcode_scores) >= 3 and Barcode_scores[-1] >= threshold and Barcode_scores[-2] >= threshold:#Check whether the sDNC fulfills robustnes criteria 85:85:85
                print 'final sDNC:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), '\n'
                print >>g, 'final sDNC:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), '\n'
                break
            else:
                npos += 1
                if npos > min([10, len(Allpos)]):
                    print 'No credible diagnosis was retrieved\n'
                    print >>g, 'No credible diagnosis was retrieved\n'
    g.close()
if __name__ == "__main__":
	main()