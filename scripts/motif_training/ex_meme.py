"""Script to execute meme"""

from subprocess import call
from itertools import chain

def run_meme(meme_path, fastafile, options, op_folder):
    """Run meme with the given options"""
    command = [meme_path + './meme']
    files = [fastafile, '-oc', op_folder]
    command += files
    others1 = ["-dna", "-mod", "zoops", "-nostatus", "-revcomp"]
    command += others1
    others2 = list(chain(*zip(options.keys(), options.values())))
    command += others2
    shell_command = ' '.join(command)
    # print(shell_command)
    call(shell_command, shell=True)
    return None

def ex_meme(meme_path, fastafile, op_folder, nmotifs=5, evt=0.001):
    options = {"-nmotifs": str(nmotifs), "-evt": str(evt)}
    run_meme(meme_path, fastafile, options, op_folder)
