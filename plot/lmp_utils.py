import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def generate_plot_filename(filename, str_title, makePDF = True):
    ''' 
    Generate a filename for a plot from the log filename with no suffix (i.e. no '.lammps')
    and a string to decribe the file
    '''
    # the file name without the .lammps or .log etc
    file_no_suf = filename.split('.')[0]
    plt_filename = file_no_suf + str_title #"_comp_plot" or "_stats_plot"
    if makePDF:
        plt_filename += ".pdf"
    else:
        plt_filename += ".png"
    return plt_filename
    
def generate_df_filename(filename):
    # the file name without the .lammps or .log etc
    file_no_suf = filename.split('.')[0]
    return file_no_suf + ".csv"
    
def get_num_structs(filename):
    num_structs = 0
    start_cnt = False
    with open(filename, 'r') as f:
        for line in f:
            if start_cnt:
                num_structs += 1
            if 'PotEng' in line:
                start_cnt = True
            if 'Loop' in line:
                start_cnt = False
    return num_structs