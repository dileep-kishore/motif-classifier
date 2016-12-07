#!/usr/bin/env python3
# @Author: Dileep Kishore <dileep>
# @Date:   2016-12-06T04:48:03-05:00
# @Filename: motif_scores.py
# @Last modified by:   dileep
# @Last modified time: 2016-12-07T01:15:44-05:00

import pandas as pd
import numpy as np
import sys

def get_bestn(fname):
    with open(fname, 'r') as fid:
        for line in fid:
            best_n = int(line.strip())
    return best_n

def get_motif_scores(fname):
    fimo_data = pd.read_table(fname)
    motif_scores = fimo_data['score']
    start, stop = fimo_data['start'], fimo_data['stop']
    start_stop = []
    for i in range(len(start)):
        start_stop.append('-'.join([str(start[i]), str(stop[i])]))
    return motif_scores, start_stop

def get_matched_seqs(fname):
    fimo_data = pd.read_table(fname)
    matched_seqs = list(fimo_data['matched sequence'])
    seq_len = len(matched_seqs[0])
    base2num = {'A':0, 'C':1, 'G':2, 'T':3}
    motif_seq_list = []
    for seq in matched_seqs:
        if seq_len != len(seq):
            print('Matched sequence all not the same length')
            sys.exit()
        else:
            motif_seq_list.append([base2num[nuc] for nuc in seq])
    motif_seq_list = np.array(motif_seq_list)
    seq_dataframe = pd.DataFrame()
    for col in range(motif_seq_list.shape[1]):
        header_name = 'seq' + str(col)
        seq_dataframe[header_name] = motif_seq_list[:,col]
    return seq_dataframe

def check_integrity(label_file, start_stop):
    label_data = pd.read_csv(label_file)
    start, stop = label_data['start'], label_data['stop']
    for i in range(len(start)):
        temp = '-'.join([str(start[i]), str(stop[i])])
        if start_stop[i] == temp:
            continue
        else:
            print('Data integrity compromised')
            sys.exit()
    return None

def extract_features(n_file, fimo_dir, labelled_file):
    best_n = get_bestn(n_file)
    fimo_file = fimo_dir + str(best_n) + '/fimo/genome/fimo.txt'
    motif_scores, start_stop = get_motif_scores(fimo_file)
    seq_dataframe = get_matched_seqs(fimo_file)
    check_integrity(labelled_file, start_stop)
    motif_info = pd.DataFrame(data={'motif_scores': motif_scores})
    return pd.concat([motif_info, seq_dataframe],axis=1)
