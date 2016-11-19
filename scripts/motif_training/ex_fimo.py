"""Script to run fimo on the entire genome/sequences"""

from subprocess import call
from itertools import chain

def run_fimo(motif_file, sequence_file, op_folder, options):
    """Run fimo with the given options"""
    command = meme_path + './fimo'
    others = list(chain(*zip(options.keys(), options.values())))
    command.append(others)
    files = ['--oc', op_folder]
    command.append(files)
    inputs = [motif_file, sequence_file]
    command.append(others)
    shell_command = ' '.join(command)
    print(shell_command)
    call(shell_command, shell=True)
    return None

def ex_fimo(motif_file, sequence_file, op_folder, thresh=0.0001<`2`>):
    options = {"--thresh": str(thresh)}
    run_fimo(motif_file, sequence_file, op_folder, options)
