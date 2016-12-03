#!/usr/bin/env python3
"""Script to extract intergenetic features from sequences"""

import pandas as pd
from numpy import isnan

def parse_intergene(intergene_file):
    intergene_data = pd.read_table(intergene_file)
    gap_start = list(intergene_data['L_END'])
    gap_stop = list(intergene_data['R_END'])
    gap_range = [(gap_start[i], gap_stop[i]) for i in range(len(gap_stop)) if not isnan(gap_start[i]) or not isnan(gap_stop[i])]
    return gap_range

def check_intergenetic(gap_range, labelled_data):
    motif_sites = pd.read_csv(labelled_data)
    motif_start = list(motif_sites['start'])
    motif_stop = list(motif_sites['stop'])
    intergenetic = [2 for _ in range(len(motif_start))]
    print(gap_range[:20])
    for i in range(len(motif_stop)):
        if any(start <= motif_start[i] <= stop for (start, stop) in gap_range):
            intergenetic[i] = 1
            continue
        if any(start <= motif_stop[i] <= stop for (start, stop) in gap_range):
            intergenetic[i] = 1
            continue
        intergenetic[i] = 0
    intergenetic_df = pd.DataFrame(data={'intergenetic': intergenetic})
    return intergenetic_df

def extract_features(intergene_file, labelled_data):
    gap_range = parse_intergene(intergene_file)
    return check_intergenetic(gap_range, labelled_data)

if __name__ == '__main__':
    fname = '../data/features/Intergenetic_Spaces.txt'
    dname = '../results/motif_training/labelled_data.csv'
    df = extract_features(fname, dname)
    df.to_csv('temp.csv', index=False)
