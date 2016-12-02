#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 01:59:51 2016

@author: Manuel Gimenez
"""

#import paths
import pandas as pd
#import numpy as np

# FIX THIS! It should be read from the package
features_data_path = "../../data/features/"

# The name of the file with the RegulonDB data of TFs and their binding site
regulondb_tfbs_filename = features_data_path + "regulondb-tfbs.txt"

# The name of the columns of the data (check the datafile for more detail)
column_names =['TF_id', 'TF_name', 'TFBS_id', 'TFBS_start', 'TFBS_end', 'DNA_strand', 'interaction_id', 'unit_regulated','effect','p_name', 'relative_center', 'TFBS_sequence', 'evidence', 'evidence_level']

# We read the data file into a Pandas Dataframe
#regulondb_df = pd.read_csv(regulondb_tfbs_filename, sep='\t', comment='#', header=None, skip_blank_lines=True)
regulondb_df = pd.read_csv(regulondb_tfbs_filename, sep='\t', comment='#', skip_blank_lines=True, names=column_names)


# Retrieves a list of TF that binds in any region between the positions defined by start_position and end_position
def get_TFs_inregion(start_position, end_position, only_strong = True):
    global regulondb_df
    df_filtered = regulondb_df[(regulondb_df['TFBS_start'] >= start_position) & (regulondb_df['TFBS_end'] <= end_position)]
    if only_strong:
        df_filtered = df_filtered[df_filtered['evidence_level'] == "Strong"]
    return df_filtered

def number_of_TFs_around_motif(motif_start, motif_end, neighbourhood_size, only_strong = True):
    region_start = motif_start - neighbourhood_size
    region_end = motif_end + neighbourhood_size
    TFs_df = get_TFs_inregion(region_start, region_end, only_strong)
    return len(TFs_df)
    
temp = get_TFs_inregion(20000,30000,False)
print(len(temp))