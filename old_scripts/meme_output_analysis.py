import os
import sys
import subprocess
import re
from collections import OrderedDict
from operator import itemgetter
from tomtom_comp import compare_motifs

class Motif:
    """Class that defines a motif"""
    def __init__(self):
        self.promoters = set()
        self.TFs = set()
        self.motif_locs = dict()
        self.evalue = dict()
        self.motif_pic = []
        self.background = ''
        self.pfm = ''

    def get_promoters(self):
        return self.promoters

    def get_TFs(self, is_eval=0):
        if is_eval:
            return self.evalue
        else:
            return self.TFs

    def add_promoters(self, new_prom, prom_loc):
        if new_prom not in self.promoters:
            self.promoters.add(new_prom)
            self.motif_locs[new_prom] = prom_loc
        return None

    def write_evalue(self, evalue, name):
        self.evalue[name] = evalue
        return None

    def add_pic(self, folname, motif_num):
        pic_name = folname + '/logo' + str(motif_num) + '.png'
        self.motif_pic.append(pic_name)
        return None

    def add_TFs(self, TF):
        self.TFs.add(TF)
        return None

    def add_background(self, background_text):
        self.background += background_text
        return None

    def add_pfm(self, pfm_text):
        self.pfm += pfm_text
        return None

    def get_background_pfm(self):
        return self.background, self.pfm

    def print_motif(self):
        print('Motif name and evalue:')
        print(self.evalue)
        print('Motif transcription factors')
        print(self.TFs)
        print('Motif promoters and their locations')
        print(self.motif_locs)
        return None

    def get_evalue(self, TF):
        return self.evalue[TF]

    def get_locs(self, promoter):
        return self.motif_locs[promoter]

    def isempty(self):
        a = not self.promoters
        b = not self.motif_locs
        c = not self.evalue
        return (a and b and c)

    # def compare_motifs(self, other):
    #     w_prom = 0
    #     comm_proms = self.promoters.intersection(other.promoters)
    #     uncomm_proms = (self.promoters | other.promoters) - comm_proms
    #     if self.promoters==other.promoters or len(uncomm_proms) < 2:
    #         if self.motif_locs == other.motif_locs:
    #             return True
    #         else:
    #             for loc in comm_proms:
    #                 x = self.motif_locs[loc]
    #                 y = other.motif_locs[loc]
    #                 #FIXME: Checking only the start site for now
    #                 # if x[1]==y[1] and -5<=x[0]-y[0]<=5 and -5<=x[2]-y[2]<=5:
    #                 #NOTE: Does strand matter?
    #                 # if -5<=x[0]-y[0]<=5 and -5<=x[2]-y[2]<=5:
    #                 if -5<=x[0]-y[0]<=5:
    #                     continue
    #                 else:
    #                     w_prom += 1
    #             if w_prom/len(comm_proms) <= 0.15:
    #                 return True
    #             else:
    #                 return False
    #     else:
    #         return False

    #FIXME: Combine background data if compared using different sets
    def combine_motifs(self, other, better_motif):
        comm_proms = self.promoters.intersection(other.promoters)
        other_proms = other.promoters - self.promoters
        #CHANGED: Average the locations of the common promoters
        for prom in comm_proms:
            for tind in [0,2]:
                self.motif_locs[prom][tind] = (self.motif_locs[prom][tind]+other.motif_locs[prom][tind])/2
        #FIXME: Code to combine promoters and locations
        for prom in other_proms:
            self.motif_locs[prom] = other.motif_locs[prom]
        if better_motif == 'target':
            self.motif_pic = self.motif_pic
        elif better_motif == 'query':
            self.motif_pic = other.motif_pic
        self.TFs = self.TFs | other.TFs
        for dict_keys in other.evalue:
            self.evalue[dict_keys] = other.evalue[dict_keys]
        return None

    def remove_TF(self):
        for prom in self.promoters:
            if 'TF' in prom:
                self.promoters.remove(prom)
                del self.motif_locs[prom]
                break
        return None

    def get_pic(self):
        picture = self.motif_pic
        if picture != []:
            return picture[len(picture)-1]
        else:
            return ' '

def get_motifs(fname, folname, text_only):
    """Identify motifs present in TFs and promoters"""
    motif_list = []
    #TODO: rewrite with while True: if fid.readline() != ''
    with open(fname, 'r') as fid:
        motif_count = 0
        curr_motif = []
        for line in fid:
            motif_start = re.search(string=line, pattern="MOTIF.*E-value = ")
            if motif_start is not None:
                # curr_motif = Motif()
                motif_list.append(curr_motif)
                evalue = line[motif_start.span()[1]:].strip()
                motif_count += 1
                if text_only:
                    tf_name = fname.split('/')[-1].split('.')[0]
                    motif_name = fname.split('/')[-1].split('.')[0] + str(motif_count)
                else:
                    tf_name = folname.split('/')[-1]
                    motif_name = tf_name + str(motif_count)
                    curr_motif.add_pic(folname, motif_count)
                curr_motif.add_TFs(tf_name)
                curr_motif.write_evalue(evalue, motif_name)
                motif_stack = []
            bck_ground = re.findall(string=line, pattern="Background letter freq")
            if len(bck_ground) == 1:
                curr_motif = Motif()
                bck_grnd_flg = 0
            if 'bck_grnd_flg' in locals() and bck_grnd_flg <= 2:
                curr_motif.add_background(line)
                bck_grnd_flg += 1
            if 'motif_stack' in locals() and line.startswith('*'):
                motif_stack.append(line)
            if 'motif_stack' in locals() and len(motif_stack) == 1:
                blk_srt = re.findall(string=line, pattern="Motif.*block diagrams")
                if len(blk_srt) == 1:
                    info_stack = []
                if 'info_stack' in locals() and line.strip() == '-'*80:
                    info_stack.append(line)
                if 'info_stack' in locals() and len(info_stack) == 1:
                    f_prom = re.findall(string=line, pattern='promoter.*_.*_')
                    f_tf = re.findall(string=line, pattern='TF.*_.*_')
                    if len(f_prom) == 1 or len(f_tf) == 1:
                        locs = line.split()[-1].split('_')
                        try:
                            locs = [int(locs[0]), locs[1], int(locs[2])]
                        except:
                            print(locs)
                        curr_motif.add_promoters(line.split()[0], locs)
                pfm_srt = re.findall(string=line, pattern=".*position-specific probability matrix")
                if len(pfm_srt) == 1:
                    pfm_stack = []
                if 'pfm_stack' in locals() and line.strip() == '-'*80:
                    pfm_stack.append(line)
                    continue
                if 'pfm_stack' in locals() and len(pfm_stack) == 1:
                    curr_motif.add_pfm(line)
    return motif_list

def write_output(out_file, motifs):
    """Write motif information into a csv file"""
    e_threshold = 0.001
    with open(out_file, 'w') as fid:
        for n_motif, motif in enumerate(motifs):
            if n_motif == 0:
                header = 'Motif,motif_pic,TFs,E-value,Promoters'
                fid.write(header)
                fid.write('\n')
            x = max(len(motif.get_promoters()), len(motif.get_TFs()))
            promoters = list(motif.get_promoters())
            TFs = list(motif.get_TFs(1).keys())
            e_list = [float(motif.get_evalue(d_TF)) for d_TF in TFs]
            e_dict = dict(zip(TFs, e_list))
            e_dict = OrderedDict(sorted(e_dict.items(), key=itemgetter(1)))
            e_items = list(e_dict.items())
            for i in range(x):
                if i < len(TFs) and e_items[i][1] > e_threshold:
                    break
                if i == 0:
                    name = 'Motif' + str(n_motif)
                    fid.write(name)
                    fid.write(',')
                    fid.write(motif.get_pic())
                    fid.write(',')
                else:
                    fid.write(',')
                    fid.write(',')
                if i < len(TFs):
                    fid.write(e_items[i][0])
                    fid.write(',')
                    fid.write(str(e_items[i][1]))
                    fid.write(',')
                else:
                    fid.write(',')
                    fid.write(',')
                if i < len(promoters):
                    fid.write(promoters[i])
                    # fid.write(',')
                    # fid.write(str(motif.get_locs(promoters[i])))
                    # fid.write(',')
                # else:
                #     fid.write(',')
                #     fid.write(',')
                fid.write('\n')

if __name__ == '__main__':
    os.chdir('../')
    if len(sys.argv) != 2:
        print('Run the program as meme_output_analysis.py <text/html>')
        sys.exit()
    text_flag = sys.argv[1]
    meme_output = 'meme_out/'
    promoter_sets = os.listdir(path=meme_output)
    promoter_sets.sort()
    out_folder = 'motif_summary/'
    try:
        subprocess.call(['rm', '-rf', out_folder])
    except:
        pass
    for prom_set in promoter_sets:
        curr_set = meme_output + prom_set + '/'
        fold_content = os.listdir(path=curr_set)
        fold_content.sort()
        motifs = []
        empty_count = 0
        for curr_cont in fold_content:
            if text_flag == 'text':
                fname = curr_set + curr_cont
                folname = ''
                motif_list = get_motifs(fname, folname, 1)
            elif text_flag == 'html':
                fname = curr_set + curr_cont + '/meme.txt'
                folname = curr_set + curr_cont
                motif_list = get_motifs(fname, folname, 0)
            good_motifs = []
            #TODO: Shuffle the motifs
            for temp in motif_list:
                if len([dum for dum in temp.get_promoters() if 'TF' in dum]) == 1:
                    temp.remove_TF()
                    good_motifs.append(temp)
            if len(motifs) == 0:
                motifs += good_motifs
                continue
            for r_motif in good_motifs:
                ind1 = 0
                is_match = 0
                r_motif.print_motif()
                if r_motif.isempty():
                    empty_count += 1
                    continue
                while ind1 < len(motifs):
                    l_motif = motifs[ind1]
                    comparision, better_motif = compare_motifs(l_motif, r_motif)
                    if comparision == True:
                        l_motif.combine_motifs(r_motif, better_motif)
                        is_match = 1
                        break
                    ind1 += 1
                if is_match == 0:
                    motifs.append(r_motif)
        print(empty_count)
        subprocess.call(['mkdir', out_folder])
        out_file = out_folder + prom_set + '.csv'
        write_output(out_file, motifs)
