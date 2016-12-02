"""Script to train a motif model from chip-seq data"""

from subprocess import call
import os
from motif_training.extract_chipdata import extract_chipdata
from motif_training.get_seqs import get_seqs
from motif_training.ex_meme import ex_meme
from motif_training.ex_fimo import ex_fimo
from motif_training.check_fitness import check_fitness
from motif_training.motif_identifier import motif_identifier
import paths
import matplotlib.pyplot as plt

def clean_up():
    results_dir = '../results/motif_training/'
    command = 'rm -rf ' + results_dir
    call(command, shell=True)
    return None

def main(TF, n, seq_len, out_dir):
    # ChIP data processing (n binding sites)
    chip_file = paths.chip_path + TF + '.xlsx'
    chip_outfile = out_dir + 'chip_out.csv'
    extract_chipdata(chip_file, n, chip_outfile)
    genome_file = paths.genome_path + 'U00096.2.fa'
    chip_fastafile = out_dir + 'chip_fasta.fa'
    get_seqs(genome_file, chip_outfile, seq_len, chip_fastafile)
    # Meme execution (n chip binding sites)
    meme_op_folder = out_dir + 'meme'
    ex_meme(paths.meme_path, chip_fastafile, meme_op_folder)
    motif_id = motif_identifier(TF, meme_op_folder)
    # Fimo execution (whole genome)
    motif_file = meme_op_folder + '/meme.txt'
    fimo_op_folder1 = out_dir + 'fimo/genome'
    ex_fimo(paths.meme_path, motif_file, genome_file, fimo_op_folder1, motif_id)
    # Fimo execution (all chip sites)
    chip_outfile_all = out_dir + 'chip_out_all.csv'
    extract_chipdata(chip_file, 'all', chip_outfile_all)
    chip_fastafile_all = out_dir + 'chip_fasta_all.fa'
    get_seqs(genome_file, chip_outfile_all, seq_len, chip_fastafile_all)
    fimo_op_folder2 = out_dir + 'fimo/chip'
    ex_fimo(paths.meme_path, motif_file, chip_fastafile_all, fimo_op_folder2, motif_id)
    # Check fitness
    if len(os.listdir(out_dir+'fimo/')) == 0:
        return None
    fit_score = check_fitness(chip_fimo=fimo_op_folder2, genome_fimo=fimo_op_folder1, chip_data=chip_outfile_all)
    return fit_score

if __name__ == '__main__':
    n = 30
    tf = 'Nac'
    seq_len = 200
    ans = input('Do you want to clean the results directory? (Y/N)\n')
    if ans == 'Y':
        clean_up()
    n_range = range(5, 100, 20)
    fitness = [0.0 for _ in n_range]
    for ind, n in enumerate(n_range):
        out_dir = paths.out_path + 'motif_training/' + str(n) + '/'
        if not os.path.exists(out_dir):
            call('mkdir -p ' + out_dir, shell=True)
        fitness[ind] = main(tf, n, seq_len, out_dir)
        try:
            print('n={0:d}, fitness={1:.3f}'.format(n, fitness[ind]))
        except TypeError:
            print('n={0:d}, fitness=None'.format(n))
    plt.plot(n_range, fitness)
    plt.show()
