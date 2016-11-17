"""Script to write a html file of the motif_summary"""

import os
import subprocess
import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader
import pdfkit

def html_writer(title, out_dir, data_set, template, rowspan_list):
    """Write the html file of the motif summary"""
    header = list(data_set.columns)
    template_vars = dict()
    template_vars['title'] = title
    template_vars['header'] = header
    template_vars['Motifs'] = data_set['Motif']
    template_vars['motif_pic'] = data_set['motif_pic']
    template_vars['TFs'] = data_set['TFs']
    template_vars['e_value'] = data_set['E-value']
    template_vars['promoters'] = data_set['Promoters']
    template_vars['row_span'] = rowspan_list
    html_op = template.render(template_vars)
    return html_op

if __name__ == '__main__':
    os.chdir('../')
    motif_dir = 'motif_summary/'
    out_dir = 'html_output/'
    env = Environment(loader=FileSystemLoader('scripts'))
    template = env.get_template('html_template.html')
    try:
        subprocess.call('rm -rf '+out_dir, shell=True)
    except:
        pass
    subprocess.call('mkdir '+out_dir, shell=True)
    motif_data = os.listdir(motif_dir)
    motif_data.sort()
    for ind, curr_set in enumerate(motif_data):
        curr_file = motif_dir + curr_set
        table_set = pd.read_csv(curr_file)
        table_set.fillna('-', inplace=1)
        title = 'set' + str(ind)
        rowspan_list = [0 for dummy in range(table_set.shape[0])]
        pic_counts = 0
        empty_counts = []
        for j, pics in enumerate(table_set['motif_pic']):
            if type(pics) == str and pics != '-':
                k = j
                cp_command = 'cp '+pics+' '+out_dir+title+'_'+str(pic_counts)+'.png'
                pic_counts += 1
                empty_counts.append(1)
                subprocess.call(cp_command, shell=True)
            if pics == '-':
                empty_counts[pic_counts-1] += 1
                rowspan_list[k] = empty_counts[pic_counts-1]
        html_data = html_writer(title, out_dir, table_set, template, rowspan_list)
        html_file = out_dir+title+'.html'
        with open(html_file, 'w') as fid:
            fid.write(html_data)
        pdf_file = html_file.split('.')[0] + '.pdf'
        pdfkit.from_file(html_file, pdf_file)
