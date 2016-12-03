#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:53:04 2016

@author: Manuel Gimenez
"""

import pandas as pd
from functools import partial

def load_from_topology_db(complete_path_to_file):
    column_names = ['wtf', 'from', 'to', 'value']
    cols_to_use = ['to', 'value']
    data_df = pd.read_csv(complete_path_to_file, sep='\t', comment='#', skip_blank_lines=True, names=column_names, usecols=cols_to_use, skiprows=1)
    return data_df

def get_DNA_topology_value(the_data, index_template, a_series):
    start = a_series[0]
    end = a_series[1]
    data_filtered = the_data[ (the_data['to'] >= start) & (the_data['to'] <= end)]
    res = data_filtered.iloc[:,1]
    labels = [index_template + str(x) for x in xrange(end-start+1) ]
    return pd.Series(res.tolist(), index=labels)

# This is the MAIN EXTRACTOR function
def get_df_topologies_of_bindingsites(positions_df):
    
    # FIX THIS! It should be read from the package paths
    features_data_path = "../../data/features/"
    
    # Filenames of the databases with DNA topology data
    mgw_db_filename = "ecoli-dnashape-mgw.db"
    helt_db_filename = "ecoli-dnashape-helt.db"
    prot_db_filename = "ecoli-dnashape-prot.db"
    roll_db_filename = "ecoli-dnashape-roll.db"
    #orchid_db_filename = "ecoli-dnashape-orchid"
    
    #We load all the DB into memory
    print("Loading DNA topology databases into memory. This can take a while...")    
    mgw_data = load_from_topology_db(features_data_path + mgw_db_filename)
    helt_data = load_from_topology_db(features_data_path + helt_db_filename)
    prot_data = load_from_topology_db(features_data_path + prot_db_filename)
    roll_data = load_from_topology_db(features_data_path + roll_db_filename)
    
    #df = pd.DataFrame({'x' : [5, 32, 50, 64], 'y' : [13, 40, 58, 72]})
    # We "currify" the function get_DNA_topology_value to use an specific db
    f_mgw = partial(get_DNA_topology_value, mgw_data, 'mgw')
    f_helt = partial(get_DNA_topology_value, helt_data, 'helt')
    f_prot = partial(get_DNA_topology_value, prot_data, 'prot')
    f_roll = partial(get_DNA_topology_value, roll_data, 'roll')
    
    # We get the topology data
    df_mgw = positions_df.apply(f_mgw, axis=1)
    df_helt = positions_df.apply(f_helt, axis=1)
    df_prot = positions_df.apply(f_prot, axis=1)
    df_roll = positions_df.apply(f_roll, axis=1)
    
    # We concatenate everything
    res = pd.concat([df_mgw, df_helt, df_prot, df_roll], axis=1)
    return res
    
# Test data frame that simualtes one with START and END positions of predicted Binding sites
#df = pd.DataFrame({'x' : [5, 32, 50, 64], 'y' : [13, 40, 58, 72]})
#temp = get_df_topologies_of_bindingsites(df)
#print(temp)