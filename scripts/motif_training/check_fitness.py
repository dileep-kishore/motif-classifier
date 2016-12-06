# @Author: Dileep Kishore <dileep>
# @Date:   2016-12-03T02:13:20-05:00
# @Filename: check_fitness.py
# @Last modified by:   dileep
# @Last modified time: 2016-12-06T14:11:06-05:00
"""Script to check the fitness/performance of a given motif"""

import pandas as pd
import numpy as np

def chip_fitness(chip_fimo, chip_file):
    """Calculate the fitness of the chip binding site predictions"""
    chip_data = pd.read_csv(chip_file)
    genes = list(chip_data['Symbol'])
    coverage = list(chip_data['Coverage'])
    hist_bins, _ = np.histogram(coverage, bins=3)
    hist_bins = list(reversed(hist_bins[]))
    weigths = []
    for ind, hist_bin in enumerate(hist_bins):
        weights += [5*(2-i)+1 for _ in range(hist_bin)]
    weights = np.array(weigths)
    norm_weigths = (weights-np.min(weights)) / (np.max(weights)-np.min(weights))
    chip_dict = dict(zip(genes, norm_weigths))
    fimo_data = pd.read_table(chip_fimo)
    num_matches = len(fimo_data)
    score = list(fimo_data['score'])
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

def check_fitness(chip_fimo, genome_fimo, chip_data):
    """Function to calculate fitness of a motif"""
    #NOTE: These parameters need to be optimized
    a = 0.7
    b = 0.3
    chip_fimo = chip_fimo + '/fimo.txt'
    genome_fimo = genome_fimo + '/fimo.txt'
    chip_matches, chip_score = chip_fitness(chip_fimo, chip_data)
    genome_matches = genome_fitness(genome_fimo)
    #TODO: Make this a weighted average wrt coverage
    norm_chip = chip_score/chip_matches
    #TODO: Only subtract intersection of genome and chip matches
    norm_match = (genome_matches-chip_matches)/chip_matches
    a = 1
    b = norm_chip
    return a*norm_chip - b*norm_match
