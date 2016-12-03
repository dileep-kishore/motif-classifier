# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 18:07:30 2016

@author: Callen
"""
'''This function returns the base frequency of the regions flanking the predicted
motif. Lm is the location of the left end of the motif. Rm is the right end of our
motif. DNA is the genomic DNA sequence. W is how far to the left and right of the motif
we want to explore.'''

def read_genome(genome_file):
    """Read genome sequence from the genome file"""
    with open(genome_file, 'r') as fid:
        seq = [line.strip() for line in fid.readlines()]
        sequence = ''.join(seq[1:])
    return sequence
def Base_Comp(Lm,Rm, DNA, W):
    #Frequency Dictionary
    Aleft=0
    Tleft=0
    Gleft=0
    Cleft=0
    Aright=0
    Tright=0
    Gright=0
    Cright=0     
    Lflank=DNA[Lm-W:Lm]
    Rflank=DNA[Rm:Rm+W]
    #I'
    for i in range(0, len(Lflank)):
        if Lflank[i]=="A":
            Aleft+=1/W
        elif Lflank[i]=="T":
            Tleft+=1/W
        elif Lflank[i]=="G":
            Gleft+=1/W
        elif Lflank[i]=="C":
            Cleft+=1/W
    for i in range(0, len(Rflank)):
        if Rflank[i]=="A":
            Aright+=1/W
        elif Rflank[i]=="T":
            Tright+=1/W
        elif Rflank[i]=="G":
            Gright+=1/W
        elif Rflank[i]=="C":
            Cright+=1/W
    return Aleft, Tleft, Gleft, Cleft, Aright, Tright, Gright, Cright

MyGenome=read_genome('U00096.2.fa.txt')
# How we can implement this function into a table/list of vectors?
# FreqMat=np.zeros((Nmotifs+1, 4))
# given a dataframe of start and end points of each motif, PandaMot, where each column represents
# frequencies of A, T, G, and C respectively.
#Also Given dataframe
# for i in range(0, Number Motifs):
#   FreqMat[i,0], FreqMat[i,1], FreqMat[i,2], FreqMat[i,3], etc= Base_Comp(PandaMot[i,0], PandaMot[i,1], DNA, W)
# Alternatively, if you want to use a list of lists, 
#try Listolist= [[0, 0, 0, 0]*Nmotifs] or a list comprehension to initialize it
