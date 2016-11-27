"""Script to train a motif model from chip-seq data"""

from motif_training.extract_chipdata import extract_chipdata
from motif_training.get_seqs import get_seqs
from motif_training.ex_meme import ex_meme
from motif_training.ex_fimo import ex_fimo
import paths

def main(TF, n, seq_len):
    # ChIP data processing (n binding sites)
    chip_file = paths.chip_path + TF + '.xlsx'
    chip_outfile = paths.out_path + 'chip_out.csv'
    extract_chipdata(chip_file, n, chip_outfile)
    genome_file = paths.genome_path + 'U00096.2.fa'
    chip_fastafile = paths.out_path + 'chip_fasta.fa'
    get_seqs(genome_file, chip_outfile, seq_len, chip_fastafile)
    # Meme execution (n chip binding sites)
    meme_op_folder = paths.out_path + 'meme'
    ex_meme(paths.meme_path, chip_fastafile, meme_op_folder)
    # Fimo execution (whole genome)
    motif_file = meme_op_folder + '/meme.txt'
    fimo_op_folder1 = paths.out_path + 'fimo/genome'
    ex_fimo(paths.meme_path, motif_file, genome_file, fimo_op_folder1)
    # Fimo execution (all chip sites)
    chip_outfile_all = paths.out_path + 'chip_out_all.csv'
    extract_chipdata(chip_file, 'all', chip_outfile_all)
    chip_fastafile_all = paths.out_path + 'chip_fasta_all.fa'
    get_seqs(genome_file, chip_outfile_all, seq_len, chip_fastafile_all)
    fimo_op_folder2 = paths.out_path + 'fimo/chip'
    ex_fimo(paths.meme_path, motif_file, chip_fastafile_all, fimo_op_folder2)
    # Check fitness
    fitness = check_fitness(chip_fimo=fimo_op_folder2, genome_fimo=fimo_op_folder1, chip_data=chip_outfile_all)
    return fitness

if __name__ == '__main__':
    n = 30
    tf = 'Nac'
    seq_len = 200
    for ind, n in enumerate(range(5, 100, 20)):
        fitness[ind] = main(tf, n, seq_len)
