#!/usr/bin/python3
"""Script to check analyzed motifs against gene expression data
    author: dileep
    date: 11-01-2016 10:30am"""

import numpy as np
import pandas as pd
import os

def read_expdata(fname, sheetnum):
    """Read expression data from excel sheet"""
    exp_data = pd.read_excel(fname, sheetname=sheetnum, skiprows=[0])
    return exp_data

def read_binddata(fname, sheetnum):
    """Reading binding data from excel sheet"""
    bind_data = pd.read_excel(fname, sheetname=sheetnum, skiprows=[0])
    return bind_data


def read_motifdata(fname):
    """Read motif data from motif summary file"""
    motif_data = pd.read_csv(fname)
    return motif_data

def parse_motifdata(motif_table):
    """Parse motif data"""
    motif_table.fillna('-', inplace=1)
    motif_dict = {}
    for ind, dat in enumerate(motif_table['Motif']):
        if type(dat) == str and dat != '-':
            curr_motif = dat
            motif_dict[curr_motif] = []
        if dat == '-':
            motif_dict[curr_motif].append(motif_table['TFs'][ind][:-1])
    return motif_dict

def parse_expdata(exp_table, tfs):
    """Parse expression data"""
    exp_dict = {'decrease': [], 'increase': [], 'no_change': []}
    not_present = []
    for tf in tfs:
        if tf not in list(exp_table['Symbol']):
            not_present.append(tf)
            continue
        fil_table = exp_table[exp_table['Symbol']==tf]
        fc = list(fil_table['Log2 FC'])[0]
        gene_id = list(fil_table['EcoCyc Locus'])[0]
        if gene_id in [d[0] for d in sum(exp_dict.values(), [])]:
            continue
        if fc < -0.05:
            exp_dict['decrease'].append([gene_id, tf, fc])
        elif fc > 0.05:
            exp_dict['increase'].append([gene_id, tf, fc])
        else:
            exp_dict['no_change'].append([gene_id, tf, fc])
    return exp_dict, not_present

if __name__ == '__main__':
    os.chdir('../')
    exp_folder = 'new_data/testbindingdata/'
    # Expression file needs to be changed for different TFs
    exp_file = exp_folder + 'Nac.xlsx'
    exp_sheet = 2
    exp_data = read_expdata(exp_file, exp_sheet)
    motif_folder = 'motif_summary/'
    # Motif file needs to be changed for each set
    motif_file = motif_folder + 'set0.csv'
    motif_data = read_motifdata(motif_file)
    motif_dict = parse_motifdata(motif_data)
    motif_exp = {}
    for motif in motif_dict:
        exp_dict, false_pos = parse_expdata(exp_data, motif_dict[motif])
        motif_exp[motif] = [exp_dict, false_pos]
