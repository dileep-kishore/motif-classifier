#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 01:59:51 2016

@author: Manuel Gimenez
"""

import pandas as pd
import numpy as np

# The name of the file with the RegulonDB data of TFs and their binding site
regulondb_tfbs_filename = "regulondb-tfbs.txt"

# The name of the columns of the data (check the datafile for more detail)
column_names =['TF_id', 'TF_name', 'TFBS_id', 'TFBS_start', 'TFBS_end', 'DNA_strand', 'interaction_id', 'unit_regulated','effect','p_name', 'relative_center', 'TFBS_sequence', 'evidence', 'evidence_level']

# We read the data file into a Pandas Dataframe
#regulondb_df = pd.read_csv(regulondb_tfbs_filename, sep='\t', comment='#', header=None, skip_blank_lines=True)
regulondb_df = pd.read_csv(regulondb_tfbs_filename, sep='\t', comment='#', skip_blank_lines=True, names=column_names)


# Retrieves a list of TF that binds in any region between the positions defined by start_position and end_position
def get_TFs_between(start_position, end_position):
    global regulondb_df
    df_filtered = regulondb_df[(regulondb_df['TFBS_start'] >= start_position) & (regulondb_df['TFBS_end'] <= end_position)]
    return df_filtered


temp = get_TFs_between(20000,30000)
print(len(temp))
    
#print(regulondb_df)