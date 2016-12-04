# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 18:07:30 2016

@author: Callen
"""
'''This function returns the base frequency of the regions flanking the predicted
motif. Lm is the location of the left end of the motif. Rm is the right end of our
motif. DNA is the genomic DNA sequence. W is how far to the left and right of the motif
we want to explore.'''

from functools import partial
import pandas as pd

def read_genome(genome_file):
    """Read genome sequence from the genome file"""
    with open(genome_file, 'r') as fid:
        seq = [line.strip() for line in fid.readlines()]
        sequence = ''.join(seq[1:])
    return sequence
def Base_Comp(Lm,Rm, DNA, W):
    
    assert (Lm - W) > 0, "DANGER: Base composition is trying to acces a negative position in the genome"
    assert (Rm + W) <= len(DNA), "DANGER: Base composition is trying to acces a outofscope position in the genome"
    
    #Frequency Dictionary
    Aleft=0.0
    Tleft=0.0
    Gleft=0.0
    Cleft=0.0
    Aright=0.0
    Tright=0.0
    Gright=0.0
    Cright=0.0   
    Lflank=DNA[Lm-W:Lm]
    Rflank=DNA[Rm:Rm+W]
    #I'
    for i in range(0, len(Lflank)):
        if Lflank[i]=="A":
            Aleft+=1.0/W
        elif Lflank[i]=="T":
            Tleft+=1.0/W
        elif Lflank[i]=="G":
            Gleft+=1.0/W
        elif Lflank[i]=="C":
            Cleft+=1.0/W
    for i in range(0, len(Rflank)):
        if Rflank[i]=="A":
            Aright+=1.0/W
        elif Rflank[i]=="T":
            Tright+=1.0/W
        elif Rflank[i]=="G":
            Gright+=1.0/W
        elif Rflank[i]=="C":
            Cright+=1.0/W
    res = [Aleft, Tleft, Gleft, Cleft, Aright, Tright, Gright, Cright]
    return res
    

def adaptor(DNA, W, a_series):
    Lm = a_series[0]
    Rm = a_series[1]
    res = Base_Comp(Lm,Rm, DNA, W)
    labels = ['compA_d', 'compT_d','compG_d','compC_d','compA_u','compT_u','compG_u','compC_u',]
    res_series = pd.Series(res, index=labels)
    return res_series
    
    
def get_df_base_composition(positions_df):
    # FIX THIS! It should be read from the package paths
    genome_data_path = "../../data/genome/"
    # Filenames of the databases with DNA topology data
    genome_filename = "U00096.2.fa"

    window_size = 25      
    genome_sequence = read_genome(genome_data_path+genome_filename)
    
    # Now we apply the function to every entry in the positions_df
    base_comp_adpated_partial = partial(adaptor, genome_sequence, window_size)
    res = positions_df.apply(base_comp_adpated_partial, axis=1)
    return res
        
df = pd.DataFrame({'x' : [50, 320, 500, 640], 'y' : [130, 400, 580, 720]})
temp = get_df_base_composition(df)
print(temp)

# How we can implement this function into a table/list of vectors?
# FreqMat=np.zeros((Nmotifs+1, 4))
# given a dataframe of start and end points of each motif, PandaMot, 
#Also Given dataframe where each column represents
# frequencies of A, T, G, and C respectively.

# for i in range(0, Number Motifs):
#   FreqMat[i,0], FreqMat[i,1], FreqMat[i,2], FreqMat[i,3], etc= Base_Comp(PandaMot[i,0], PandaMot[i,1], DNA, W)
# Alternatively, if you want to use a list of lists, 
#try Listolist= [[0, 0, 0, 0]*Nmotifs] or a list comprehension to initialize it
