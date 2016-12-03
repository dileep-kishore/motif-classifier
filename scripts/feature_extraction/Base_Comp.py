# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 18:07:30 2016

@author: Callen
"""
'''This function returns the base frequency of the regions flanking the predicted
motif. Lm is the location of the left end of the motif. Rm is the right end of our
motif. DNA is the genomic DNA sequence. W is how far to the left and right of the motif
we want to explore.'''
def Base_Comp(Lm,Rm, DNA, W):
    #Frequency Dictionary
    A=0
    T=0
    G=0
    C=0    
    Lflank=DNA[Lm-W:Lm]
    Rflank=DNA[Rm:Rm+W]
    for i in range(0, len(Lflank)):
        if DNA[Lflank[i]]=="A":
            A+=1/W
        elif DNA[Lflank[i]]=="T":
            T+=1/W
        elif DNA[Lflank[i]]=="G":
            G+=1/W
        elif DNA[Lflank[i]]=="C":
            C+=1/W
    return A,T,G,C