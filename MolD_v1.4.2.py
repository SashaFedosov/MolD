"""
This script compiles rDNC-based DNA diagnoses for a pre-defined taxa in a dataset. This is the MAIN WORKING VERSION v1.4.2 including latest bug fixes as of July 2024.
This version already implements the new functionalities:
- automatic trimming of the alignment to match the median sequence length (should seq len distribution and provided NumberN settings require that)
- expanded qCLADEs setting
- diagnoses in pairwise comparisons of taxa.
- selection of a reference sequence for indexing DNCs

THIS version LACKS SPART compatibility

"""
import os, sys
import random
import argparse
from io import StringIO
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

#***STEP 4 ANALYSIS OF OUTPUT rDNCs
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

#SPECIFIC FUNCTIONS FOR THE rDNCs
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
    #print(["ROUND", len(SEQ)*Pdiff/100, round(len(SEQ)*Pdiff/100)])
    n = round(len(SEQ)*Pdiff/100)
    N = random.sample(list(range(1, n)), 1)[0]
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
            iSTS = random.sample(list(range(len(seqlist))), nseq_per_clade_to_screw)
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

#NEWFUNCYIONS
def medianSeqLen(listofseqs):#OCT2022
    seqlens = [i.count('A')+i.count('C')+i.count('G')+i.count('T') for i in listofseqs]
    medlen = sorted(seqlens)[int(len(seqlens)/2)]
    medseq = listofseqs[seqlens.index(medlen)]
    start = min([medseq.find('A'),medseq.find('C'),medseq.find('G'),medseq.count('T')])
    if not 'N' in medseq[start:]:
        end = len(medseq)
    else:
        for i in range(start, len(medseq), 1):
            if medseq[i:].count('N') == len(medseq[i:]):
                end = i
                break
    return medlen, start, end

def getAllPairs(taxalist):
    uniquetaxapairs = []
    stl = sorted(list(set(taxalist)))
    for i in range(len(stl)):
        for j in range(i+1, len(stl)):
            uniquetaxapairs.append([stl[i],stl[j]])
    return uniquetaxapairs
        

################################################READ IN PARAMETER FILE AND DATA FILE

def get_args(): #arguments needed to give to this script
    parser = argparse.ArgumentParser(description="run MolD")
    required = parser.add_argument_group("required arguments")
    required.add_argument("-i", help="textfile with parameters of the analysis", required=True)
    return parser.parse_args()

def mainprocessing(gapsaschars=None, taxalist=None, taxonrank=None, cutoff=None, numnucl=None, numiter=None, maxlenraw=None, maxlenrefined=None, iref=None, pdiff=None, nmax=None, thresh=None, tmpfname=None, origfname=None):
    ParDict = {}
    if not(all([gapsaschars, taxalist, taxonrank, cutoff, numnucl, numiter, maxlenraw, maxlenrefined, iref, pdiff, nmax, thresh, tmpfname])):
        args = get_args()
        with open(args.i) as params:
            for line in params:
                line = line.strip()
                if line.startswith('#'):
                    pass
                else:
                    if len(line.split('=')) == 2 and len(line.split('=')[1]) != 0:
                        ParDict[line.split('=')[0]] = line.split('=')[1].replace(' ', '')#VERYNEW
    else:
        ParDict['Gaps_as_chars'] = gapsaschars
        ParDict['qTAXA'] = taxalist
        ParDict['Taxon_rank'] = taxonrank
        ParDict['INPUT_FILE'] = tmpfname
        ParDict['ORIG_FNAME'] = origfname# I do not understand this ORIG_FNAME, it throws an error with the command line
        ParDict['Cutoff'] = cutoff
        ParDict['NumberN'] = numnucl
        ParDict['Number_of_iterations'] = numiter
        ParDict['MaxLen1'] = maxlenraw
        ParDict['MaxLen2'] = maxlenrefined
        ParDict['Iref'] = iref
        ParDict['Pdiff'] = pdiff
        #ParDict['PrSeq'] = prseq
        ParDict['NMaxSeq'] = nmax
        ParDict['Scoring'] = thresh
        ParDict['OUTPUT_FILE'] = "str"        
    print(ParDict)
    ############################################# #VERYNEW HOW GAPS ARE TREATED
    #REQUIRES A NEW FIELD IN THE GUI
    if ParDict['Gaps_as_chars'] == 'yes':
        gaps2D = True#VERYNEW
    else:#VERYNEW
        gaps2D = False#VERYNEW
    
    ############################################ DATA FILE
    checkread = open(ParDict['INPUT_FILE'], 'r')
    firstline = checkread.readline()
    checkread.close()
    imported=[]#set up a new list with species and identifiers
    if '>' in firstline:
        f = open(ParDict['INPUT_FILE'], 'r') 
        for line in f:#VERYNEW - THE DATA READING FROM ALIGNMENT FILE IS ALL REVISED UNTIL f.close()
            line=line.rstrip()
            if line.startswith('>'):
                data = line[1:].split('|')
                if len(data) != 2:
                    print('Check number of entries in', data[0])
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
    else:
        f = open(ParDict['INPUT_FILE'], 'r') 
        for line in f:#VERYNEW - THE DATA READING FROM ALIGNMENT FILE IS ALL REVISED UNTIL f.close()
            line=line.rstrip()
            data = line.split('\t')
            if len(data) != 3:
                print('Check number of entries in', data[0])
            ID = data[0]
            thetax = data[1]
            if gaps2D == True:#
                DNA = data[2].upper().replace('-', 'D').replace('R', 'N').replace('Y', 'N').replace('S','N').replace('W','N').replace('M','N').replace('K','N')#VERYNEW
            else:
                DNA = data[2].upper().replace('-', 'N').replace('R', 'N').replace('Y', 'N').replace('S','N').replace('W','N').replace('M','N').replace('K','N')#VERYNEW
            imported.append([ID, thetax, DNA])
        f.close()
    
    if 'NumberN' in list(ParDict.keys()):#How many ambiguously called nucleotides are allowed
        NumberN = int(ParDict['NumberN'])
    else:
        NumberN = 5
    
    if len(set([len(i[2]) for i in imported])) != 1:
        print('Alignment contains sequences of different lengths:', set([len(i[2]) for i in imported]))
    else:
        mlen, sstart, send = medianSeqLen([i[2] for i in imported])#OCT2022 - start
        if mlen + NumberN < len(imported[0][2]):
            Slice = True
            FragmentLen = mlen
            corr = sstart
        else:
            Slice = False
            FragmentLen = len(imported[0][2])
            sstart = 0
            send = FragmentLen+1
            corr = 0
    raw_records=[]
    for i in imported:
        if ParDict['Iref'] != 'NO' and ParDict['Iref'].split(',')[0] == i[0] and ParDict['Iref'].split(',')[1] in ['ex', 'excl', 'out']:
            continue
        else:
            if i[2][sstart:send].count('N') < NumberN and len(i[2][sstart:send]) == FragmentLen:
                newi = [i[0], i[1], i[2][sstart:send]]
                raw_records.append(newi)#OCT2022 - end
    print('\n########################## PARAMETERS ######################\n')#VERYNEW
    #print('input file:', ParDict['ORIG_FNAME']) #Outcommented ORIG_FNAME
    print('input file:', ParDict['INPUT_FILE']) #Replacement of the line above
    print('Coding gaps as characters:', gaps2D)
    print('Maximum undetermined nucleotides allowed:', NumberN)
    print('Length of the alignment:', len(imported[0][2]),'->', FragmentLen)
    print('Indexing reference:', ParDict['Iref'].replace('NO', 'Not set').replace('in', 'included').replace('ex', 'excluded'))
    print('Read in', len(raw_records), 'sequences')
    
    PosArrays, VarPosList = PositionArrays([i[2] for i in raw_records])#VERYNEW
   
    #############################################READ IN OTHER ANALYSIS PARAMETERS
    ##OCT2022 - start
    withplus = []
    P2 = []
    shift = True
    if ParDict['qTAXA'][0] == '>':#THIS OPTION DIAGNOSES ALL TAXA WITH MORE THAN USER-DEFINED NUMBER OF SEQUENCES AVAILABLE
        NumSeq = int(ParDict['qTAXA'][1:])
        Taxarecords = [i[1] for i in raw_records]
        qCLADEs = []
        for j in set(Taxarecords):
            if Taxarecords.count(j) >= NumSeq:
                qCLADEs.append(j)
    elif ParDict['qTAXA'].startswith('P:'):#THIS OPTION DIAGNOSES ALL TAXA CONTAINING A USER-DEFINED PATTERN IN THE NAME
        pattern = ParDict['qTAXA'].split(':')[1]
        Taxarecords = [i[1] for i in raw_records]
        qCLADEs = []
        for j in set(Taxarecords):
            if pattern in j:
                qCLADEs.append(j)
    elif ParDict['qTAXA'].startswith('P+:'):#THIS OPTION POOLS ALL TAXA CONTAINING A USER-DEFINED PATTERN IN THE NAME IN ONE TAXON AND DIAGNOSES IT
        pattern = ParDict['qTAXA'].split(':')[1]
        Taxarecords = set([i[1] for i in raw_records if pattern in i[1]])
        spp = '+'.join(sorted(list(Taxarecords)))
        qCLADEs = [spp]
        nrecords = []
        for rec in raw_records:
            if rec[1] in Taxarecords:
                nrecords.append([rec[0], spp, rec[2]])
            else:
                nrecords.append(rec)
        raw_records = nrecords
    else:#THIS OPTION DIAGNOSES ALL TAXA FROM A USER-DEFINED LIST; TAXA MAY BE COMBINED BY USING '+'
        qCLADEs = []
        allrecs = ParDict['qTAXA'].split(',')
        for item in allrecs:
            if item in ['ALL', 'All', 'all']:#THIS OPTION DIAGNOSES ALL TAXA IN THE DATASET
                qCLADEs = list(set([i[1] for i in raw_records]))
            elif item in [i[1] for i in raw_records]:
                qCLADEs.append(item)
            elif '+' in item:
                withplus.append(item)
            elif 'VS' in item:
                P2.append(item.split('VS'))
            else:
                print('UNRECOGNIZED TAXON', item)
    #OCT2022 - end
    print('query taxa:', len(qCLADEs+withplus), '-', str(sorted(qCLADEs)+sorted(withplus)).replace('[','').replace(']','').replace("'", ''))#1.3
    
    if 'Cutoff' in list(ParDict.keys()):#CUTOFF Number of the informative positions to be considered, default 100
        Cutoff = ParDict['Cutoff']#VERYNEW
    else:
        Cutoff = 100
    print('Cutoff set as:', Cutoff)
    if 'Number_of_iterations' in list(ParDict.keys()):#Number iterations of MolD
        N1 = int(ParDict['Number_of_iterations'])
    else:
        N1 = 10000
    print('Number iterations of MolD set as:', N1)
    
    if 'MaxLen1' in list(ParDict.keys()):#Maximum length for the raw mDNCs
        MaxLen1 = int(ParDict['MaxLen1'])
    else:
        MaxLen1 = 12
    print('Maximum length of raw mDNCs set as:', MaxLen1)
    
    if 'MaxLen2' in list(ParDict.keys()):#Maximum length for the refined mDNCs
        MaxLen2 = int(ParDict['MaxLen2'])
    else:
        MaxLen2 = 7
    print('Maximum length of refined mDNCs set as:', MaxLen2)
    
    if 'Pdiff' in list(ParDict.keys()):#Percent difference
        Percent_difference = float(ParDict['Pdiff'])
    else:
        if int(ParDict['Taxon_rank']) == 1:#read in taxon rank to configure Pdiff parameter of artificial dataset
            Percent_difference = 1
        else:
            Percent_difference = 5
    print('simulated sequences up to', Percent_difference, 'percent divergent from original ones')
    
    if 'NMaxSeq' in list(ParDict.keys()):#Maximum number of sequences per taxon to be modified
        Seq_per_clade_to_screw = int(ParDict['NMaxSeq'])
    else:
        Seq_per_clade_to_screw = 10####!changed value
    print('Maximum number of sequences modified per clade', Seq_per_clade_to_screw)
    
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
    #print(ParDict['Scoring'], 'scoring of the rDNCs; threshold in two consequtive runs:', threshold)
    print('scoring of the rDNCs; threshold in two consequtive runs:', threshold)
    #OCT2022 - start
    if corr > 1 and len(ParDict['Iref'].split(',')) == 2:
        print('\nNOTE: The alignment was trimmed automatically to match median sequences length. The analysed slice starts from the site',str(sstart+1),'and ends on the site',str(send+1),'. The site indexing in the DNCs as in the provided reference.')
    if corr > 1 and ParDict['Iref'] == 'NO':
        print('\nNOTE: The alignment was trimmed automatically to match median sequences length. The analysed slice starts from the site',str(sstart+1),'and ends on the site',str(send+1),'. The site indexing in the rDNC as in the untrimmed alignment.')
    #OCT2022 - end
    thephrase = 'The DNA diagnosis for the taxon'
    ###################################################IMPLEMENTATION
    #Setting up a new class just for the convenient output formatting
    class SortedDisplayDict(dict):#this is only to get a likable formatting of the barcode
        def __str__(self):
            return "[" + ", ".join("%r: %r" % (key, self[key]) for key in sorted(self)) + "]"
    
    class SortedDisplayDictVerbose(dict):#this is only to get a likable formatting of the barcode
        def __str__(self):
            return ", ".join("%r %r" % (self[key],'in the site '+str(key)) for key in sorted(self)).replace("'", '').replace("A", "'A'").replace("C", "'C'").replace("G", "'G'").replace("T", "'T'")+'.'
    
    #Calling functions and outputing results
    if ParDict['OUTPUT_FILE'] == "str":
        g = StringIO()
    else:
        g = open(ParDict['OUTPUT_FILE'], "w")#Initiating output file
    #VERYNEW
    print('<h4>########################## PARAMETERS ######################</h4>', file=g)
    #print("<p>", 'input file:', ParDict['ORIG_FNAME'], "</p>", file=g) #outcommented ORIG_FNAME
    print("<p>", 'input file:', ParDict['INPUT_FILE'], "</p>", file=g)
    print("<p>", 'Coding gaps as characters:', gaps2D, "</p>", file=g)
    print("<p>", 'Maximum undetermined nucleotides allowed:', NumberN, "</p>", file=g)
    print("<p>", 'Length of the alignment:', len(imported[0][2]),'->', FragmentLen, "</p>", file=g)
    print("<p>", 'Indexing reference:', ParDict['Iref'].replace('NO', 'Not set').replace('in', 'included').replace('ex', 'excluded'), "</p>", file=g)#OCT2022
    print("<p>", 'Read in', len(raw_records), 'sequences', "</p>", file=g)
    print("<p>", 'query taxa:', len(qCLADEs+withplus), '-', str(sorted(qCLADEs)+sorted(withplus)).replace('[','').replace(']','').replace("'", ''), "</p>", file=g)#1.3
    print("<p>", 'Cutoff set as:', Cutoff, "</p>", file=g)
    print("<p>", 'Number iterations of MolD set as:', N1, "</p>", file=g)
    print("<p>", 'Maximum length of raw mDNCs set as:', MaxLen1, "</p>", file=g)
    print("<p>", 'Maximum length of refined mDNCs set as:', MaxLen2, "</p>", file=g)
    print("<p>", 'simulated sequences up to', Percent_difference, 'percent divergent from original ones', "</p>", file=g)
    print("<p>", 'Maximum number of sequences modified per clade', Seq_per_clade_to_screw, "</p>", file=g)
    #print("<p>", ParDict['Scoring'], 'scoring of the rDNCs; threshold in two consequtive runs:', threshold, "</p>", file=g)
    print("<p>", 'scoring of the rDNCs; threshold in two consequtive runs:', threshold, "</p>", file=g)
    if corr > 1 and len(ParDict['Iref'].split(',')) == 2:
        print('<h4>NOTE: The alignment was trimmed automatically to match median sequences length. The analysed slice starts from the site',str(sstart+1),'and ends on the site',str(send+1),'. The site indexing in the DNCs as in the provided reference.</h4>', file=g)
    if corr > 1 and ParDict['Iref'] == 'NO':
        print('<h4>NOTE: The alignment was trimmed automatically to match median sequences length. The analysed slice starts from the site',str(sstart+1),'and ends on the site',str(send+1),'. The site indexing in the rDNC as in the untrimmed alignment.</h4>', file=g)

    print('<h4>########################### RESULTS ##########################</h4>', file=g)
    
    for qCLADE in sorted(qCLADEs) + sorted(withplus):#OCT2022
        if '+' in qCLADE:
            if shift == True:
                old_records = raw_records
                shift == False
            else:
                raw_records = old_records
            spp = qCLADE.split('+')
            nrecords = []
            for rec in raw_records:
                if rec[1] in spp:
                    nrecords.append([rec[0], qCLADE, rec[2]])
                else:
                    nrecords.append(rec)
            raw_records = nrecords
        print('\n**************', qCLADE, '**************')
        print('<h4>**************', qCLADE, '**************</h4>', file=g)
        Clades, clade_sorted_seqs, shared_positions = Step1(raw_records)#STEP1
        x,y,z,pures = C_VP_PP(clade_sorted_seqs, qCLADE, shared_positions, Cutoff)#STEP2 ####! added pures
        newy = {key:y[key] for key in y if not key in pures} ####! newline
        print('Sequences analyzed:', len(clade_sorted_seqs[qCLADE]))
        print("<p>",'Sequences analyzed:', len(clade_sorted_seqs[qCLADE]), "</p>", file=g)
        ND_combinations = [[item] for item in pures] ####! before ND_combinations were initiated as an empty list
        print('single nucleotide mDNCs:', len(pures), '-', str(SortedDisplayDict({pos+corr : y[pos-1] for pos in [i+1 for i in pures]}))[1:-1])#OCT2022
        print("<p>",'single nucleotide mDNCs*:',len(pures), '-', str(SortedDisplayDict({pos+corr : y[pos-1] for pos in [i+1 for i in pures]}))[1:-1], "</p>", file=g)#OCT2022
        N = 1 ####!
        while N > 0:#STEP3
            try:
                q = Diagnostic_combinations(qCLADE, x, newy, N1, MaxLen1, MaxLen2) ####! newy instead of y
            except IndexError:
                print(N, 'IndexError')
                continue
            for comb in q:
                if not comb in ND_combinations:
                    ND_combinations.append(comb)
            N-=1
        ND_combinations.sort(key=len)
        #################################### mDNC output
        try:
            Nind, KeyPos = IndependentKey(ND_combinations)#STEP4
        except IndexError:
            print('no mDNCs recovered for', qCLADE)#VERYNEW
            print("<p>", 'no mDNCs recovered for', "</p>", qCLADE, file=g)#VERYNEW
            continue
        Allpos = []#Create list of all positions involved in mDNCs
        for comb in ND_combinations:
            for pos in comb:
                if not pos in Allpos:
                    Allpos.append(pos)
        print('\nmDNCs retrieved:', str(len(ND_combinations)) + '; Sites involved:', str(len(Allpos)) + '; Independent mDNCs:', len(Nind))#VERYNEW
        print("<p>", 'mDNCs* retrieved:', str(len(ND_combinations)) + '; Sites involved:', str(len(Allpos)) + '; Independent mDNCs**:', len(Nind), "</p>", file=g)#VERYNEW
        print('Shortest retrieved mDNC:', SortedDisplayDict({pos+corr : y[pos-1] for pos in [i+1 for i in ND_combinations[0]]}), '\n')#OCT2022
        print("<p>",'Shortest retrieved mDNC*:', SortedDisplayDict({pos+corr : y[pos-1] for pos in [i+1 for i in ND_combinations[0]]}), "</p>", file=g)#OCT2022
        ######################################################## rDNC output
        Barcode_scores = []#Initiate a list for rDNC scores
        npos = len(ND_combinations[0])
        BestBarcode = 'none'####! newline
        while npos <= min([10, len(Allpos)]):#in this loop the positions are added one-by-one to a rDNC and the rDNC is then rated on the artificially generated datasets
            Barcode = GenerateBarcode_new(ND_combinations, npos)#Initiate a rDNC
            Barcode_score = 0#Initiate a score to rate a rDNC
            N = 100
            while N > 0:
                NComplist, NCPP = Screwed_dataset_new(raw_records, Seq_per_clade_to_screw, PosArrays, VarPosList, Percent_difference, qCLADE, Cutoff)#Create an artificial dataset VERYNEW
                NBarcode = [i for i in Barcode if i in list(NCPP.keys())]
                if len(Barcode) - len(NBarcode) <= 1 and ConditionD(NBarcode, NComplist, NCPP) == True:####! new condition (first) added
                    Barcode_score +=1
                N -=1
            print(npos, 'rDNC_score (100):', [k+corr+1 for k in Barcode], '-', Barcode_score)#VERYNEW
            print("<p>", npos, 'rDNC_score (100):', [k+corr+1 for k in Barcode], '-', Barcode_score, "</p>", file=g)#OCT2022
            if Barcode_score >= threshold and len(Barcode_scores) == 1: ###1.3
                BestBarcode = Barcode###1.3
            if Barcode_score >= threshold and len(Barcode_scores) > 1 and Barcode_score >= max(Barcode_scores): ###1.3
                BestBarcode = Barcode####!newline
            Barcode_scores.append(Barcode_score)
            if len(Barcode_scores) >= 2 and Barcode_scores[-1] >= threshold and Barcode_scores[-2] >= threshold:#Check whether the rDNC fulfills robustnes criteria 85:85:85
                print('Final rDNC:', SortedDisplayDict({pos+corr : y[pos-1] for pos in [i+1 for i in BestBarcode]}))#OCT2022
                print("<p>",'Final rDNC***:', SortedDisplayDict({pos+corr : y[pos-1] for pos in [i+1 for i in BestBarcode]}), "</p>", file=g)#OCT2022
                print('\n',thephrase, qCLADE, 'is:', SortedDisplayDictVerbose({pos+corr : y[pos-1] for pos in [i+1 for i in BestBarcode]}))#OCT2022
                print("<p>", thephrase, qCLADE, 'is:', SortedDisplayDictVerbose({pos+corr : y[pos-1] for pos in [i+1 for i in BestBarcode]}), "</p>", file=g)#OCT2022
                break
            else:# VERY NEW FROM HERE ONWARDS
                npos += 1
                if npos > min([10, len(Allpos)]):
                    if BestBarcode != 'none':
                        print('The highest scoring rDNC for taxon', qCLADE, 'is:', SortedDisplayDictVerbose({pos+corr : y[pos-1] for pos in [i+1 for i in BestBarcode]}))#OCT2022
                        print("<p>", 'The highest scoring rDNC*** for taxon', qCLADE, 'is:', SortedDisplayDictVerbose({pos+corr : y[pos-1] for pos in [i+1 for i in BestBarcode]}), "</p>", file=g)#OCT2022
                    else:
                        print('No sufficiently robust DNA diagnosis for taxon', qCLADE, 'was retrieved')
                        print("<p>", 'No sufficiently robust DNA diagnosis for taxon', qCLADE, 'was retrieved', "</p>", file=g)
    #OCT2022 - start
    print("<h4>", '################################# EXPLANATIONS ####################################', "</h4>", file=g)
    print("<p>", '  * mDNC -(=minimal Diagnostic nucleotide combination) is a combination of nucleotides at specified sites of the alignment,', "</p>", file=g)
    print("<p>", '    unique for a query taxon. Therefore it is sufficient to differentiate a query taxon from all reference taxa in a dataset.', "</p>", file=g)
    print("<p>", '    Because it comprises minimal necessary number of nucleotide sites to differentiate a query, any mutation in the mDNC in'  "</p>", file=g)
    print("<p>", '    single specimen of a query taxon will automatically disqualify it as a diagnostic combination.', "</p>", file=g)
    print("<p>",  "</p>", file=g)
    print("<p>", ' ** two or more mDNCs are INDEPENDENT if they constitute non-overlapping sets of nucleotide sites.', "</p>", file=g)
    print("<p>",  "</p>", file=g)
    print("<p>", '*** rDNC -(=robust/redundant Diagnostic nucleotide combination) is a combination of nucleotides at specified sites of the alignment,', "</p>", file=g)
    print("<p>", '    unique for a query taxon and (likewise mDNC) sufficient to differentiate a query taxon from all reference taxa in a dataset.', "</p>", file=g)
    print("<p>", '    However, rDNC comprises more than a minimal necessary number of diagnostic sites, and therefore is robust to single nucleotide', "</p>", file=g)
    print("<p>", '    replacements. Even if a mutation arises in one of the rDNC sites, the remaining ones will (with high probability) remain sufficient ', "</p>", file=g)
    print("<p>", '    to diagnose the query taxon', "</p>", file=g)
    print("<h4>", '    Final diagnosis corresponds to rDNC', "</h4>", file=g)
    #OCT2022 - end
    
    if ParDict['OUTPUT_FILE'] == "str":
        contents = g.getvalue()
        os.unlink(ParDict['INPUT_FILE'])
    else:
        contents = None
    g.close()
    #OCT2022 - start
    if len(P2) != 0:
        ext = '.'+ParDict['OUTPUT_FILE'].split('.')[-1]
        h = open(ParDict['OUTPUT_FILE'].replace(ext, '_pairwise'+ext), "w")#Initiating output file
        if len(withplus) != 0:
            raw_records = old_records
        taxain = [i[1] for i in raw_records]
        tpairs = []
        for alist in P2:
            if alist.count('all')+alist.count('All')+alist.count('ALL') == 2:
                tpairs = getAllPairs(taxain)
            elif alist.count('all')+alist.count('All')+alist.count('ALL') == 1:
                thetax = [i for i in taxain if i in alist][0]
                print(thetax)
                for atax in sorted(list(set(taxain))):
                    if atax != thetax:
                        tpairs.append([thetax, atax])
            else:
                for apair in getAllPairs(alist):
                    tpairs.append(apair)
        for apair in tpairs:
            t1 = apair[0]
            t2 = apair[1]
            p2records = [i for i in raw_records if i[1] in [t1, t2]]
            print('\n**************', t1, 'VS', t2,'**************')
            print('<h4>**************', t1, 'VS', t2, '**************</h4>', file=h)
            C2, css2, sp2 = Step1(p2records)#STEP1
            x2,y2,z2,pures2 = C_VP_PP(css2, t1, sp2, '>0')#STEP2 ####! added pures
            Pairphrase = 'Each of the following '+ str(len(pures2))+' sites is invariant across sequences of '+ t1+ ' and differentiates it from '+ t2+': '
            counterPures = {}
            for site in pures2:
                counterPures[site] = "'or'".join(list(set([thing[site] for thing in css2[t2] if thing[site] != 'N'])))
                Pairphrase = Pairphrase + str(site+corr+1)+" ('"+str(y2[site])+"' vs '"+str(counterPures[site])+"'), "
            print(Pairphrase[:-2].replace("''","'-'"))
            print("<p>",Pairphrase[:-2].replace("''","'-'"),'</h4>', file=h)#OCT2022
            x2r,y2r,z2r,pures2r = C_VP_PP(css2, t2, sp2, '>0')#STEP2 ####! added pures
            Pairphraser = 'Each of the following '+ str(len(pures2r))+' sites is invariant across sequences of '+ t2+ ' and differentiates it from '+ t1+': '
            counterPuresr = {}
            for site in pures2r:
                counterPuresr[site] = "'or'".join(list(set([thing[site] for thing in css2[t1] if thing[site] != 'N'])))
                Pairphraser = Pairphraser + str(site+corr+1)+" ('"+str(y2r[site])+"' vs '"+str(counterPuresr[site])+"'), "
            print(Pairphraser[:-2].replace("''","'-'"))
            print("<p>",Pairphraser[:-2].replace("''","'-'"),'</h4>', file=h)#OCT2022
        h.close()
            
    #OCT2022 - end
    return contents, qCLADEs
if __name__ == "__main__":
    c, q = mainprocessing()
  
