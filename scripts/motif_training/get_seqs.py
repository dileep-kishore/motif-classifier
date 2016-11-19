"""Script to get nucleotide sequences from chip-seq positions"""

def read_genome(genome_file):
    """Read genome sequence from the genome file"""
    with open(genome_file, 'r') as fid:
        seq = [line.strip() for line in fid.readlines()]
        sequence = ''.join(seq[1:])
    return sequence

def extract_seqs(genome, positions, seq_len):
    """Get chip sequences from the genome"""
    seq_list = []
    for pos in positions:
        start = pos - seq_len//2
        end = pos + seq_len//2
        seq_list.append(genome[start:end])
    return seq_list

def write_seqs(chip_data, seqs, fastafile):
    """Write sequences into a fasta file"""
    with open(fastafile, 'w') as fid:
        for ind, gene in enumerate(list(chip_data['Symbol'])):
            header = '>'+gene+'\n'
            fid.write(header)
            fid.write(seqs[ind])
            fid.write('\n')
    return None

def get_seqs(genome_file, chip_data, seq_len, fastafile):
    """Extract chip sequences from the genome"""
    genome = read_genome(genome_file)
    positions = list(chip_data['Position'])
    sequences = extract_seqs(genome, positions, seq_len)
    write_seqs(chip_data, sequences, fastafile)
