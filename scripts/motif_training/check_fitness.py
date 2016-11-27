"""Script to check the fitness/performance of a given motif"""

import pandas as pd

def chip_fitness(chip_fimo, chip_file):
    """Calculate the fitness of the chip binding site predictions"""
    fimo_data = pd.read_table(chip_fimo)
    chip_data = pd.read_csv(chip_file)
    num_matches = len(fimo_data)
    score = list(fimo_data['score'])
    seq_name = list(fimo_data['sequence name'])
    fimo_dict = dict(zip(seq_name, score))
    genes = list(chip_data['Symbol'])
    tot_score = 0
    for seq in fimo_dict:
        tot_score += fimo_dict[seq]
    return num_matches, tot_score

def genome_fitness(genome_fimo):
    """Calculate the fitness of the genome binding site predictions"""
    fimo_data = pd.read_table(genom_fimo)
    fitness = len(fimo_data)
    return fitness

def check_fitness(chip_fimo, genome_fimo, chip_data):
    """Function to calculate fitness of a motif"""
    a = 
    b = 
    chip_mat, chip_score = chip_fitness(chip_fimo, chip_data)
    genome_matches = genome_fitness(genome_fimo)
    return a*(chip_score/chip_matches) - b*(genome_matches-chip_matches)
