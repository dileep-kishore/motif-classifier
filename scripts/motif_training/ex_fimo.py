# @Author: Dileep Kishore <dileep>
# @Date:   2016-12-01T20:47:32-05:00
# @Filename: ex_fimo.py
# @Last modified by:   dileep
# @Last modified time: 2016-12-06T22:56:04-05:00



"""Script to run fimo on the entire genome/sequences"""

from subprocess import call
from itertools import chain

def run_fimo(meme_path, motif_file, sequence_file, op_folder, options):
    """Run fimo with the given options"""
    command = [meme_path + './fimo']
    others = list(chain(*zip(options.keys(), options.values())))
    command += others
    files = ['--oc', op_folder]
    command += files
    inputs = [motif_file, sequence_file]
    command += inputs
    shell_command = ' '.join(command)
    # print(shell_command)
    call(shell_command, shell=True)
    return None

def ex_fimo(meme_path, motif_file, sequence_file, op_folder, motif_id, thresh=0.0001):
    options = {"--thresh": str(thresh), "--verbosity": str(1), "--motif": motif_id}
    fimo_dir = '/'.join(op_folder.split('/')[:-1])
    try:
        call('mkdir '+fimo_dir, shell=True)
    except:
        pass
    run_fimo(meme_path, motif_file, sequence_file, op_folder, options)
