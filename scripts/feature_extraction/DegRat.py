# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 18:25:24 2016

@author: Callen
"""

'''This function gives the ratio of the degrees of W base pairs/turn in the 
flanking regions to the degrees of W base pairs, given that BDNA is usually
about 10.55 base pairs per turn. Lm is left end. Rm is right end. W is 
how far beyond the motif we want to explore. DNA is our genomic sequence'''

def DegRat(Lm, Rm, DNA, W, DegData):
    LTotDeg=0.
    RTotDeg=0.
    Lflank=DNA[Lm-W:Lm]
    Rflank=DNA[Rm:Rm+W]
    for i in range(0, len(Lflank)):
        LTotDeg+=DegData[Lflank[i]]
        RTotDeg+=DegData[Rflank[i]]
    RDegRat=RTotDeg/(360./10.55 * len(Rflank))
    LDegRat=RTotDeg/(360./10.55 * len(Lflank))
    return RDegRat, LDegRat