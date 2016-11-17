"""Script to generate output from meme"""
# Note: Uses only these two plugins in addition to meme (might run on server)
import subprocess
import os
import sys

def main(input_folder, output_folder, text_only):
    """Identify common motifs in TF-promoters and promoter seqs using meme"""
    # part = ['forward', 'reverse']
    fol_in = input_folder + '/'
    fol_out = output_folder + '/'
    subprocess.call(["mkdir", fol_out])
    in_files = os.listdir(fol_in)
    in_files.sort()
    for count, fname in enumerate(in_files):
        print(str(count)+'/'+str(len(in_files)), end='\r')
        in_file = fol_in + fname
        out_folder = fol_out + fname[:-3]
        if text_only:
            out_file = out_folder + '.txt'
            shell_command = ["/home/dileep/meme/bin/./meme", in_file, "-oc", out_folder, "-dna", "-mod", "zoops", "-nmotifs", "5", "-nostatus", "-revcomp", "-evt", "0.001", "-text", ">", out_file]
        else:
            shell_command = ["/home/dileep/meme/bin/./meme", in_file, "-o", out_folder, "-dna", "-mod", "zoops", "-nmotifs", "5", "-nostatus", "-revcomp", "-evt", "0.001"]
        subprocess.call(' '.join(shell_command), shell=True)
    return None

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Run the program as get_meme_output.py <text/html>')
        sys.exit()
    text_flag = sys.argv[1]
    os.chdir('../')
    fol_in = 'meme_in/'
    fol_out = 'meme_out/'
    answer = input("Do you want to delete the output folder Y/N?")
    if answer == 'Y':
        try:
            subprocess.call(["rm", "-rf", fol_out])
            subprocess.call(["mkdir", fol_out])
        except:
            pass
    promoter_sets = os.listdir(fol_in)
    promoter_sets.sort()
    subprocess.call(["mkdir", fol_out])
    for prom_set in promoter_sets:
        print("Analyzing: ", prom_set)
        if text_flag == 'text':
            main(fol_in+prom_set, fol_out+prom_set, 1)
        elif text_flag == 'html':
            main(fol_in+prom_set, fol_out+prom_set, 0)
