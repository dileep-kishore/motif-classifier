# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 18:07:30 2016

@author: Callen
"""
'''This function returns the base frequency of the regions flanking the predicted
motif. Lm is the location of the left end of the motif. Rm is the right end of our
motif. DNA is the genomic DNA sequence. W is how far to the left and right of the motif
we want to explore.'''
def Base_Comp(Lm,Rm, DNA, W, DNA):
    #Frequency Dictionary
    Freqy={"A":0, "T":0, "G":0,  "C":0,}
    Lflank=DNA[Lm-W:Lm]
    Rflank=DNA[Rm:Rm+W]
    for i in range(0, len(Lflank)):
        if DNA[Lflank[i]]=="A":
            Freqy["A"]+=1/W
        elif DNA[Lflank[i]]=="T":
            Freqy["T"]+=1/W
        elif DNA[Lflank[i]]=="G":
            Freqy["G"]+=1/W
        elif DNA[Lflank[i]]=="C":
            Freqy["C"]+=1/W
    return Freqy