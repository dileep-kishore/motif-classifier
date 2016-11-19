"""Script to train a motif model from chip-seq data"""

from motif_training.extract_chipdata import extract_chipdata
from motif_training.get_seqs import get_seqs
from motif_training.ex_meme import ex_meme
from motif_training.ex_fimo import ex_fimo
from paths import path_init

def main(TF, n, seq_len):
    chip_file = chip_path + TF + '.xlsx'
    chip_outfile = out_path + 'chip_out.csv'
    extract_chipdata(chip_file, n, chip_outfile)
    genome_file = genome_path + 'U00096.2.fa'
    chip_fastafile = out_path + 'chip_fasta.fa'
    get_seqs(genome_file, chip_outfile, seq_len, chip_fastafile)
    meme_op_folder = out_path + 'meme'
    ex_meme(chip_fastafile, meme_op_folder)
    motif_file = meme_op_folder + '/meme.txt'
    fimo_op_folder = out_path + 'fimo'
    ex_fimo(motif_file, genome_file, fimo_op_folder)

if __name__ == '__main__':
    path_init()
    main('Nac',n,seq_len)
