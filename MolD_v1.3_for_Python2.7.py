"""
This script compiles sDNC-based DNA diagnoses for a pre-defined taxa in a dataset. This is the MAIN WORKING VERSION

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
    pures = []####! newline
    for key in CShN:
        newcomplist = []
        for k in complist:
            if k[key] == CShN[key]:
                newcomplist.append(k)
            else: continue
        cutoffs[key] = len(complist) - len(newcomplist)
        if len(newcomplist) == 0:####! newline
            pures.append(key)####! newline
    CPP = []
    for key in sorted(cutoffs, key = cutoffs.get, reverse = True):
        CPP.append(key)
    if CUTOFF[0] == '>':#VERYNEW
        Clade_priority_positions = {pos:CShN[pos] for pos in CPP if cutoffs[pos] > int(CUTOFF[1:])}#VERYNEW
    else:#VERYNEW
        Clade_priority_positions = {}
        for position in CPP[:int(CUTOFF)]:#Here you define how many of the clade shared combinations are used in subsequent search
            Clade_priority_positions[position] = CShN[position]
    return complist, Clade_priority_positions, cutoffs, pures####! pures added

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
            if j[pos] == CPP[pos] or j[pos] == 'N':#VERYNEW
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
        raw_comb = step_reduction_complist(qCLADE, complist, CPP, Achecked_ind)
        if len(raw_comb) <= maxlen1:
            refined_comb = RemoveRedundantPositions(raw_comb, complist, CPP)
            if len(refined_comb) <= maxlen2 and sorted(refined_comb) not in bestlists:
                bestlists.append(sorted(refined_comb))
        n=n-1
    bestlists.sort(key=len)
    return bestlists

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
def PositionArrays(Motifs):#VERYNEW ALL FUNCTION NEW
    PositionArrays = []
    VarPosList = []
    for i in range(len(Motifs[0])):
        Const = True
        array = [Motifs[0][i]]
        for j in range(len(Motifs[1:])):
            if Motifs[j][i] != 'N':
                if Motifs[j][i] != array[-1]:
                    Const = False
                array.append(Motifs[j][i])
        PositionArrays.append(array)
        if Const == False:
            VarPosList.append(i)
    return PositionArrays, VarPosList

def random_sequence_new(SEQ, PositionArrays, VarPosList, Pdiff):#VERYNEW FUNCTION REVISED
    n = len(SEQ)*Pdiff/100
    N = random.sample(range(1, n), 1)[0]
    PosToChange = random.sample([p for p in VarPosList if SEQ[p] != 'D'], N)#this is a new definition to keep alignment gaps unchanged
    NEWSEQ = ''
    for i in range(len(SEQ)):
        if i not in PosToChange:
            NEWSEQ = NEWSEQ + SEQ[i]
        else:
            newarray = [j for j in PositionArrays[i] if j != SEQ[i]]
            newbase = random.sample(newarray, 1)[0]
            NEWSEQ = NEWSEQ + newbase
    return NEWSEQ

def GenerateBarcode_new(Diagnostic_combinations, length):#VERYNEW FUNCTION REVISED - This function calculates diagnostic combinations and assembles a barcode of desired length for a query taxon. First all single position DNCs are added, then based on the frequency of a nucleotide position in the DNCs of the 2 positions, and then based on the frequency of a position in longer DNCs
    len1 = []
    len2 = []
    lenmore = []
    for comb in Diagnostic_combinations:
        if len(comb) == len(Diagnostic_combinations[0]):
            for i in comb:
                len1.append(i)
        elif len(comb) == len(Diagnostic_combinations[0])+1:
            for j in comb:
                len2.append(j)
        else:
            for k in comb:
                lenmore.append(k)
    if len(Diagnostic_combinations[0]) == 1:
        Setin = len1
    else:
        Setin = []
        for pos in sorted(len1, key=len1.count, reverse = True):
            if not pos in Setin:
                Setin.append(pos)
    for pos1 in sorted(len2, key=len2.count, reverse = True):
        if not pos1 in Setin:
            Setin.append(pos1)
    for pos2 in sorted(lenmore, key=lenmore.count, reverse = True):
        if not pos2 in Setin:
            Setin.append(pos2)
    return Setin[:length]

def Screwed_dataset_new(raw_records, nseq_per_clade_to_screw, PositionArrays, VarPosList, Percent_difference, Taxon, Cutoff):#VERYNEW FUNCTION REVISED
    clades=[]
    for i in range(len(raw_records)):
        Clade=raw_records[i][1]
        if Clade not in clades:
            clades.append(Clade)
    clade_sorted_seqs = {}
    for letter in clades:
        clade_sorted_seqs[letter]=[]
        for i in range(len(raw_records)):
            if raw_records[i][1]==letter:
                clade_sorted_seqs[letter].append(raw_records[i][2])
    
    for clade in clades:
        seqlist = clade_sorted_seqs[clade]
        newseqs = []
        if len(seqlist) > nseq_per_clade_to_screw:
            iSTS = random.sample(range(len(seqlist)), nseq_per_clade_to_screw)
            for k in range(len(seqlist)):
                if k in iSTS:
                    newseq = random_sequence_new(seqlist[k], PositionArrays, VarPosList, Percent_difference)
                else:
                    newseq = seqlist[k]
                newseqs.append(newseq)
        elif len(clade_sorted_seqs[clade]) == nseq_per_clade_to_screw:
            for k in range(len(seqlist)):
                newseq = random_sequence_new(seqlist[k], PositionArrays, VarPosList, Percent_difference)
                newseqs.append(newseq)
        else:
            for i in range(nseq_per_clade_to_screw):
                seq = random.sample(seqlist, 1)[0]
                newseq = random_sequence_new(seq, PositionArrays, VarPosList, Percent_difference)
                newseqs.append(newseq)
        clade_sorted_seqs[clade] = newseqs
    
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
                    ParDict[line.split('=')[0]] = line.split('=')[1].replace(' ', '')#VERYNEW

    ############################################# #VERYNEW HOW GAPS ARE TREATED
    #REQUIRES A NEW FIELD IN THE GUI
    if ParDict['Gaps_as_chars'] == 'yes':
        gaps2D = True#VERYNEW
    else:#VERYNEW
        gaps2D = False#VERYNEW
    
    ############################################ DATA FILE    
    f = open(ParDict['INPUT_FILE'], 'r') 
    imported=[]#set up a new list with species and identifiers
    for line in f:#VERYNEW - THE DATA READING FROM ALIGNMENT FILE IS ALL REVISED UNTIL f.close()
        line=line.rstrip()
        if line.startswith('>'):
            data = line[1:].split('|')
            if len(data) != 2:
                print 'Check number of entries in', data[0]
                #break
            data.append('')
            imported.append(data)
        else:
            if gaps2D == True:#
                DNA = line.upper().replace('-', 'D').replace('R', 'N').replace('Y', 'N').replace('S','N').replace('W','N').replace('M','N').replace('K','N')#VERYNEW
            else:
                DNA = line.upper().replace('-', 'N').replace('R', 'N').replace('Y', 'N').replace('S','N').replace('W','N').replace('M','N').replace('K','N')#VERYNEW
            imported[-1][-1] = imported[-1][-1]+DNA
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
    print '\n########################## PARAMETERS ######################\n'#VERYNEW
    print 'input file:', ParDict['INPUT_FILE'] #VERYNEW
    print 'Coding gaps as characters:', gaps2D
    print 'Maximum undetermined nucleotides allowed:', NumberN
    print 'Length of the alignment:', FragmentLen    
    print 'Read in', len(raw_records), 'sequences'
    
    PosArrays, VarPosList = PositionArrays([i[2] for i in raw_records])#VERYNEW
   
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
    print 'query taxa:', len(qCLADEs), '-', str(sorted(qCLADEs)).replace('[','').replace(']','').replace("'", '')#1.3
    
    if 'Cutoff' in list(ParDict.keys()):#CUTOFF Number of the informative positions to be considered, default 100
        Cutoff = ParDict['Cutoff']#VERYNEW
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
    
    if 'Pdiff' in list(ParDict.keys()):#Percent difference
        Percent_difference = int(ParDict['Pdiff'])
    else:
        if int(ParDict['Taxon_rank']) == 1:#read in taxon rank to configure Pdiff parameter of artificial dataset
            Percent_difference = 2
        else:
            Percent_difference = 5
    print 'simulated sequences up to', Percent_difference, 'percent divergent from original ones'
    
    
    if 'NMaxSeq' in list(ParDict.keys()):#Maximum number of sequences per taxon to be modified
        Seq_per_clade_to_screw = int(ParDict['NMaxSeq'])
    else:
        Seq_per_clade_to_screw = 10####!changed value
    print 'Maximum number of sequences modified per clade', Seq_per_clade_to_screw
    
    if 'Scoring' in list(ParDict.keys()):
        if ParDict['Scoring'] == 'lousy':
            threshold = 66####!changed value
        elif ParDict['Scoring'] == 'moderate':
            threshold = 75####!changed value
        elif ParDict['Scoring'] == 'stringent':
            threshold = 90####!changed value
        elif ParDict['Scoring'] == 'very_stringent':
            threshold = 95####!changed value
        else:
            threshold = 75####!changed value
    else:
        threshold = 75####!changed value
    print ParDict['Scoring'], 'scoring of the sDNCs; threshold in two consequtive runs:', threshold
    
    ###################################################IMPLEMENTATION
    #Setting up a new class just for the convenient output formatting
    class SortedDisplayDict(dict):#this is only to get a likable formatting of the barcode
        def __str__(self):
            return "[" + ", ".join("%r: %r" % (key, self[key]) for key in sorted(self)) + "]"
    
    #Calling functions and outputing results
    g = open(ParDict['OUTPUT_FILE'], "w")#Initiating output file
    #VERYNEW
    print >>g, '########################## PARAMETERS ######################'
    print >>g, 'input file:', ParDict['INPUT_FILE'] #VERYNEW
    print >>g, 'Coding gaps as characters:', gaps2D
    print >>g, 'Maximum undetermined nucleotides allowed:', NumberN
    print >>g, 'Length of the alignment:', FragmentLen    
    print >>g, 'Read in', len(raw_records), 'sequences'
    print >>g, 'query taxa:', len(qCLADEs), '-', str(sorted(qCLADEs)).replace('[','').replace(']','').replace("'", '')#1.3
    print >>g, 'Cutoff set as:', Cutoff
    print >>g, 'Number iterations of MolD set as:', N1
    print >>g, 'Maximum length of raw pDNCs set as:', MaxLen1
    print >>g, 'Maximum length of refined pDNCs set as:', MaxLen2
    print >>g, 'simulated sequences up to', Percent_difference, 'percent divergent from original ones'
    print >>g, 'Maximum number of sequences modified per clade', Seq_per_clade_to_screw
    print >>g, ParDict['Scoring'], 'scoring of the sDNCs; threshold in two consequtive runs:', threshold
    print >>g, '\n########################### RESULTS ##########################'
    
    for qCLADE in sorted(qCLADEs):#1.3
        print '\n**************', qCLADE, '**************'
        print >>g, '\n**************', qCLADE, '**************'
        Clades, clade_sorted_seqs, shared_positions = Step1(raw_records)#STEP1
        x,y,z,pures = C_VP_PP(clade_sorted_seqs, qCLADE, shared_positions, Cutoff)#STEP2 ####! added pures
        newy = {key:y[key] for key in y if not key in pures} ####! newline
        print 'Sequences analyzed:', len(clade_sorted_seqs[qCLADE])
        print >>g, 'Sequences analyzed:', len(clade_sorted_seqs[qCLADE])
        ND_combinations = [[item] for item in pures] ####! before ND_combinations were initiated as an empty list
        print 'single nucleotide pDNCs:', len(pures), '-', str(SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in pures]}))[1:-1]#VERYNEW
        print >>g, 'single nucleotide pDNCs:',len(pures), '-', str(SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in pures]}))[1:-1]#VERYNEW
        N = 1 ####!
        while N > 0:#STEP3
            try:
                q = Diagnostic_combinations(qCLADE, x, newy, N1, MaxLen1, MaxLen2) ####! newy instead of y
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
            print 'no pDNCs recovered for', qCLADE#VERYNEW
            print >>g, 'no pDNCs recovered for', qCLADE#VERYNEW
            continue
        Allpos = []#Create list of all positions involved in pDNCs
        for comb in ND_combinations:
            for pos in comb:
                if not pos in Allpos:
                    Allpos.append(pos)
        print 'pDNCs retrieved :', str(len(ND_combinations)) + '; Positions involved:', str(len(Allpos)) + '; Independent pDNCs:', len(Nind)#VERYNEW
        print >>g, 'pDNCs retrieved :', str(len(ND_combinations)) + '; Positions involved:', str(len(Allpos)) + '; Independent pDNCs:', len(Nind)#VERYNEW
        print 'Shortest retrieved diagnostic combination:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in ND_combinations[0]]})
        print >>g, 'Shortest retrieved diagnostic combination:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in ND_combinations[0]]})
        ######################################################## sDNC output
        Barcode_scores = []#Initiate a list for sDNC scores
        npos = len(ND_combinations[0])
        BestBarcode = 'none'####! newline
        while npos <= min([10, len(Allpos)]):#in this loop the positions are added one-by-one to a sDNC and the sDNC is then rated on the artificially generated datasets
            Barcode = GenerateBarcode_new(ND_combinations, npos)#Initiate a sDNC
            Barcode_score = 0#Initiate a score to rate a sDNC
            N = 100
            while N > 0:
                NComplist, NCPP = Screwed_dataset_new(raw_records, Seq_per_clade_to_screw, PosArrays, VarPosList, Percent_difference, qCLADE, Cutoff)#Create an artificial dataset VERYNEW
                NBarcode = [i for i in Barcode if i in list(NCPP.keys())]
                if len(Barcode) - len(NBarcode) <= 1 and ConditionD(NBarcode, NComplist, NCPP) == True:####! new condition (first) added
                    Barcode_score +=1
                N -=1
            print npos, 'sDNC_score (100):', [k+1 for k in Barcode], '-', Barcode_score#VERYNEW
            print >>g, npos, 'sDNC_score (100):', [k+1 for k in Barcode], '-', Barcode_score#VERYNEW
            if Barcode_score >= threshold and len(Barcode_scores) == 1: ###1.3
                BestBarcode = Barcode###1.3
            if Barcode_score >= threshold and len(Barcode_scores) > 1 and Barcode_score >= max(Barcode_scores): ####! newline
                BestBarcode = Barcode####!newline
            Barcode_scores.append(Barcode_score)
            if len(Barcode_scores) >= 2 and Barcode_scores[-1] >= threshold and Barcode_scores[-2] >= threshold:#Check whether the sDNC fulfills robustnes criteria 85:85:85
                print 'final sDNC:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in BestBarcode]})
                print >>g, 'final sDNC:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in BestBarcode]})
                break
            else:# VERY NEW FROM HERE ONWARDS
                npos += 1
                if npos > min([10, len(Allpos)]):
                    if BestBarcode != 'none':
                        print 'The highest scoring sDNC:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in BestBarcode]})####!newline
                        print >>g, 'The highest scoring sDNC:', SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in BestBarcode]})####!newline
                    else:
                        print 'No sufficiently robust diagnosis was retrieved'
                        print >>g, 'No sufficiently robust was retrieved'
    g.close()
if __name__ == "__main__":
	main()
