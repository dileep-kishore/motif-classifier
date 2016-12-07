# @Author: Dileep Kishore <dileep>
# @Date:   2016-12-03T02:13:20-05:00
# @Filename: check_fitness.py
# @Last modified by:   dileep
# @Last modified time: 2016-12-07T00:31:14-05:00
"""Script to check the fitness/performance of a given motif"""

import pandas as pd
import numpy as np
from Bio import motifs

def chip_fitness(chip_fimo, chip_file):
    """Calculate the fitness of the chip binding site predictions"""
    chip_data = pd.read_csv(chip_file)
    genes = list(chip_data['Symbol'])
    coverage = list(chip_data['Coverage'])
    hist_bins, _ = np.histogram(coverage, bins=3)
    hist_bins = list(reversed(hist_bins))
    weights = []
    for ind, hist_bin in enumerate(hist_bins):
        weights += [5*(2-ind)+1 for _ in range(hist_bin)]
    weights = np.array(weights)
    norm_weigths = (weights) / (np.max(weights))
    chip_dict = dict(zip(genes, norm_weigths))
    fimo_data = pd.read_table(chip_fimo)
    num_matches = len(fimo_data)
    scores = list(fimo_data['score'])
    seq_name = list(fimo_data['sequence name'])
    fimo_dict = dict(zip(seq_name, scores))
    tot_score = 0
    for seq in fimo_dict:
        tot_score += fimo_dict[seq] * chip_dict[seq]
    return num_matches, tot_score

def genome_fitness(genome_fimo):
    """Calculate the fitness of the genome binding site predictions"""
    fimo_data = pd.read_table(genome_fimo)
    fitness = len(fimo_data)
    return fitness

def get_chipseq_ranges(chip_fimo, chip_all, motif_file):
    chip_all_data = pd.read_csv(chip_all)
    chip_fimo_data = pd.read_table(chip_fimo)
    motif_id = set(chip_fimo_data['#pattern name']).pop()
    with open(motif_file, 'r') as fid:
        motif_record = motifs.parse(fid, 'MEME')
    mot_len = int(motif_record[motif_id-1].length)
    mot_range = []
    for pos in list(chip_all_data['Position']):
        mot_range.append((pos-mot_len-10, pos+mot_len+10))
    return mot_range

def label_genome_data(genome_fimo, motif_range):
    genome_fimo_data = pd.read_table(genome_fimo)
    mot_start = list(genome_fimo_data["start"])
    mot_stop = list(genome_fimo_data["stop"])
    label = ['' for _ in range(len(mot_stop))]
    for i in range(len(mot_start)):
        if any(start <= mot_start[i] <= stop for (start, stop) in motif_range):
            label[i] = 1
            continue
        if any(start <= mot_stop[i] <= stop for (start, stop) in motif_range):
            label[i] = 1
            continue
        label[i] = 0
    return label

def check_fitness(chip_fimo, genome_fimo, chip_data, motif_file):
    """Function to calculate fitness of a motif"""
    #NOTE: These parameters need to be optimized
    a = 3
    b = 1
    chip_fimo = chip_fimo + '/fimo.txt'
    genome_fimo = genome_fimo + '/fimo.txt'
    chip_matches, chip_score = chip_fitness(chip_fimo, chip_data)
    genome_matches = genome_fitness(genome_fimo)
    #CHANGED: This is now a weighted average wrt coverage
    # norm_chip = chip_score/chip_matches
    norm_chip = chip_score
    #TODO: Only subtract intersection of genome and chip matches
    motif_range = get_chipseq_ranges(chip_fimo, chip_data, motif_file)
    labels = label_genome_data(genome_fimo, motif_range)
    norm_match = sum([1 for d in labels if d == 0])
    print('normchip={0:.3f}, normmatch={1:.3f}'.format(norm_chip, norm_match))
    return a*norm_chip - b*norm_match
