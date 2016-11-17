#!/usr/bin/python3
"""Script to compare motifs using the tomtom module in meme
    author: dileep
    date: 10-25-2016 5:50am"""

import numpy as np
import pandas as pd
from subprocess import call
# from meme_output_analysis import Motif

def write_tomtom_input(fname, background, pfm_list):
    """Write motif file in tomtom input format"""
    with open(fname, 'w') as fid:
        fid.write('MEME version 4\n\n')
        fid.write('ALPHABET= ACGT\n\n')
        fid.write('strands: + -\n\n')
        fid.write(background)
        fid.write('\n')
        for pfm in pfm_list:
            fid.write(pfm)
            fid.write('\n')

def get_pfm(r_motif_list, alphabet):
    """Get position-specific probability matrix"""
    alpha_ind = dict(zip(alphabet, range(0, len(alphabet))))
    pfm = np.zeros((len(r_motif_list[0]),len(alphabet)), dtype=float)
    for s_ind in range(len(r_motif_list)):
        for m_ind in range(len(r_motif_list[s_ind])):
            nuc = r_motif_list[s_ind][m_ind]
            pfm[m_ind, alpha_ind[nuc]] += 1/len(r_motif_list)
    return pfm

def generate_random_motifs(background, pfm):
    """Generate random motifs to make up the target database"""
    background = background.split('\n')[1].split(' ')
    #NOTE: prob are string
    bck_dict = dict(zip(background[::2], background[1::2]))
    alphabet = list(bck_dict.keys())
    prob_list = []
    for key in alphabet:
        prob_list.append(float(bck_dict[key]))
    width = int(pfm.split('\n')[1].split(' ')[5])
    r_motif_list = []
    for _ in range(20):
        r_motif = np.random.choice(alphabet, size=width, replace=True, p=prob_list)
        r_motif_list.append(r_motif)
    r_pfm = get_pfm(r_motif_list, alphabet)
    pfm_preface = pfm.split('\n')[1] + '\n'
    pfm_text = ''
    for row in range(r_pfm.shape[0]):
        for col in range(r_pfm.shape[1]):
            temp = " {:.6f} ".format(r_pfm[row, col])
            pfm_text += temp
        pfm_text += '\n'
    return pfm_preface + pfm_text

def parse_tomtom_output(fname):
    tom_table = pd.read_table(fname)
    target_table = tom_table[tom_table['Target ID']=='target']
    return target_table['E-value'], target_table['Overlap'], target_table['Optimal offset']

def compare_motifs(query_motif, target_motif):
    """Compare query motif and targer motif"""
    temp_folder = 'tomtom/'
    call('rm -rf '+temp_folder, shell=True)
    call('mkdir '+temp_folder, shell=True)
    query_file = temp_folder + 'query.txt'
    target_file = temp_folder + 'target.txt'
    tom_input = temp_folder + 'tom_input.txt'
    target_dbsize = 50
    q_pfm_list = []
    # Query motif
    q_bck, q_pfm = query_motif.get_background_pfm()
    q_pfm = 'MOTIF query\n' + q_pfm
    q_pfm_list.append(q_pfm)
    write_tomtom_input(query_file, q_bck, q_pfm_list)
    # Target motif
    t_pfm_list = []
    t_bck, t_pfm = target_motif.get_background_pfm()
    # print(q_pfm)
    t_pfm = 'MOTIF target\n' + t_pfm
    t_pfm_list.append(t_pfm)
    for t_ind in range(1, target_dbsize):
        r_pfm = generate_random_motifs(t_bck, t_pfm)
        r_pfm = 'MOTIF random' + str(t_ind) +'\n' + r_pfm
        t_pfm_list.append(r_pfm)
    write_tomtom_input(target_file, t_bck, t_pfm_list)
    tomtom = '~/meme/bin/tomtom '
    inputs = query_file + ' ' + target_file
    output = ' -o ' + temp_folder + 'tomout/'
    verbosity = ' -verbosity 1'
    tomtom += inputs + output + verbosity
    call(tomtom, shell=True)
    output_file = temp_folder + 'tomout/tomtom.txt'
    evalue, overlap, offset = parse_tomtom_output(output_file)
    evalue = list(evalue)[0]
    overlap = list(overlap)[0]
    offset = list(offset)[0]
    per_overlap = overlap / (len(q_pfm.split('\n'))-3)
    # print(per_overlap, evalue, offset)
    if per_overlap > 0.75 and evalue < 0.0001:
        if offset <= 0 and per_overlap > 0.99:
            return True, 'target'
        elif offset > 0 and per_overlap > 0.99:
            return True, 'query'
        else:
            return True, 'target'
    else:
        return False, 'none'
