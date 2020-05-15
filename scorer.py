
### FUNCTION: GET PFM
### Gets the PFM for a TF binding site from the .pfm file

def GetPFM(PFM_file):
    MyFile = open(PFM_file, 'r')
    code = ['A', 'C', 'G', 'T']
    freq = []
    for line in MyFile:
        if not line.startswith('>'):
            # Take the line for the nucleotide and split it in a list
            # containing the freq of that nucleotide for each position
            tmp = line.strip().split(' ')
            
            while '' in tmp:
                tmp.remove('') #cleans the list of empty spaces
             # Converts each string into a float object (decimal number)
             # and then appends the whole thing to the freq list 
                
            freq.append([float(x) for x in tmp])
    MyFile.close()
    PFM = dict(zip(code, freq)) #creates a dictionary matching nt and freqs
    return PFM

### FUNCTION: PFM to PSSM
### Calculates the PSSM for that PFM, including pseudocounts
    
def PFMtoPSSM(PFM, pseudo = 1):
    
    import math
    # We calculate the total absolute frequency for each position.
    # Here we create a list containing the sum of each position over all
    # possible nucleotides. PFM.values() contains a list of four lists, each
    # for every nucleotide, and the function zip pairs up all the first
    # elements of that lists, then the second ones... into a zip object
    # that can be used as an iterable by the function sum(x). The *
    # is added so *list(PFM.values()) elements can be taken as arguments 
    # for the zip function
    
    total_freq = [sum(x) for x in zip(*list(PFM.values()))]
    PSSM = {} #creates an empty dictionary
    
    # Now, for each one of the nucleotides, the PSSM value is calculated
    # First, we loop through every nt in PFM.keys(). We add an entry to the
    # dictionary PSSM containing a key being that nucleotide. 
    # Then, the PSSM value is associated to that key.
    
    # The PSSM value is calculated pairing the absolute frequency of that
    # nt in that position with the total absolute frequency of that position
    # throught the zip function. Then, the log2 logarightm of the relative
    # frequency (absolute frequency + pseudocount / total absolute frequency)
    # is divided by 0.25 (the probability of finding that nucleotide by chance)
    # Basically, this is a log(odds ratio)
    
    for nt in PFM.keys():
        PSSM[nt] = [math.log2(((x+pseudo)/y)/0.25) 
                    for x, y in zip(PFM[nt], total_freq)]
    
    PSSM['N'] = [min(x) for x in zip(*list(PSSM.values()))]
    
    #adds a 'N' row, containing the minimum values (so relative score is 0)
    
    return PSSM
### FUNCTION: REVCOMP
### Returns the reverse complement of a nucleotide sequence
    
def RevComp(Seq):
    comp = str.maketrans('ACGT', 'TGCA') #make a translation table
    return Seq[::-1].translate(comp) #reverse the sequence and give its comp

### FUNCTION: GETSEQ
### Gets a nucleotidic sequence from a FASTA file and returns it uppercase

def GetSeq(FASTA_file):
    MyFile = open(FASTA_file, 'r')
    Seq = ""
    for line in MyFile:
        if not line.startswith('>'):
            Seq += line.strip()
    MyFile.close()
    return(Seq.upper())

### FUNCTION: GETSCORE
### Taken a Seq and a PSSM, calculates the score of each nucleotide
### of that sequence being a match for the PSSM
    
def GetScore(Seq, PSSM):
    
    Seq = Seq.upper()
    ScoreDict = {'+':{}, '-':{}}
    
    for sequence in [Seq, RevComp(Seq)]: #gets the strandness
        
        if sequence == Seq:
            strand = '+'
        else:
            strand = '-'
        
        for nt in range(len(sequence)): #foreach nucleotide of the sequence
            
            subseq = sequence[nt:(nt+len(PSSM['A']))]
            
            if len(subseq) == len(PSSM['A']): #checks it has not reach the end 
                
                site_score = 0
                    
                for pos in range(len(subseq)): #foreach nt of the subseq
                   
                    # Get the nt of that position and then the score of that
                    # nucleotide for that position 
                    site_score += PSSM[subseq[pos]][pos]
                    
                ScoreDict[strand][str(nt)] = site_score 
            
            # we keep nt 1 being nt 0 because 
            # its convinient for localization later 
    
    return(ScoreDict)

### FUNCTION: GETRELATIVESCORE
### Given a Seq and a PSSM, calculates the relative score of each nucleotide
### being the start of a match for the PSSM
    
def GetRelativeScore(Seq, PSSM):
    
    # First, we get max and min possible scores for a sequence in a PSSM
    # Basically, we sum up all the max values and min values for each position,
    # getting the scores for the best sequence and the worst one
    
    MaxScore = sum([max(x) for x in zip(*list(PSSM.values()))])
    MinScore = sum([min(x) for x in zip(*list(PSSM.values()))])
    
    ScoreDict = GetScore(Seq, PSSM) #get a score Dict, sorted by strand
    RelScoreList = {'+':{}, '-':{}} #creates an identical one for the RelScore
    
    for strand in ScoreDict:
        for nt in ScoreDict[strand].keys():
            # for each nucleotide, get the relative score above 0
            # we "substract" the MinScore because its negative (- x - = +),
            # so we are adding that "baseline" to get every score as a positive
            # number 
            
            RelScoreList[strand][nt] = (
                (ScoreDict[strand][nt] - MinScore) / (MaxScore - MinScore)
                )
    
    return RelScoreList

            

if __name__ == '__main__':
   
    import sys
    
    MyPFM = GetPFM("MA0150.1.pfm") # this can be changed for each PFM
    print(MyPFM)
    MyPSSM = PFMtoPSSM(MyPFM)
    print(MyPSSM)

    MySeq = GetSeq(sys.argv[1]) #FASTA_file seq as first argument. FULL PATH
    print(MySeq)
    MyScore = GetScore(MySeq, MyPSSM)
    MyRelScore = GetRelativeScore(MySeq, MyPSSM)
    
    outputfile = 'output/'+sys.argv[1].split("\\")[-1]+"_output"
    # gets the filename splitting the fullpath and conserving the last element
    MyFile = open(outputfile, 'w')
    Output = "Position\tScore\tRelative Score\tStrand\tMotif\tSite\n"
    
    for i in MyScore.keys():
        for f, g in zip(MyScore[i].items(), MyRelScore[i].items()):
            Output += (f[0] + #Position
                       "\t" + str(round(f[1], 3)) + #Score
                       "\t" + str(round(g[1], 3)) + #Relative Score
                       "\t")
            
            if i == '+':
                motif = (MySeq[
                    (int(f[0])):(int(f[0]) + len(MyPSSM['A']))
                    ])
            if i == '-':
                motif = (RevComp(MySeq)[
                    (int(f[0])):(int(f[0]) + len(MyPSSM['A']))
                    ])
            
            Output += (i + '\t' 
                       + motif + '\t' + sys.argv[1].split("\\")[-1] + "\n")
    MyFile.write(Output)
    MyFile.close()
