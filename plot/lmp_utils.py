import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def generate_plot_filename(filename, str_title, makePDF = True):
    ''' 
    Generate a filename for a plot from the log filename with no suffix (i.e. no '.lammps')
    and a string to decribe the file
    '''
    # the file name without the .lammps or .log etc
    # but there might be decimals so only cut out the last one
    file_no_suf = '.'.join( filename.split('.')[:-1] )
    plt_filename = file_no_suf + str_title #"_comp_plot" or "_stats_plot"
    if makePDF:
        plt_filename += ".pdf"
    else:
        plt_filename += ".png"
    return plt_filename
    
def generate_df_filename(filename):
    # the file name without the .lammps or .log etc
    file_no_suf = '.'.join( filename.split('.')[:-1] )
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
    
def convert_to_ase(filename):
    '''
        the full xyz filename
        Converts the types (1,2,3,4) to (Si, C, Si, C)
        to be properly read by ASE's xyz reader
    '''
    file_no_suf = '.'.join( filename.split('.')[:-1] )
    file_ase = file_no_suf + "_ase.xyz"
    open(file_ase,"w").close()
    outfile = open(file_ase,"w")
     
    start = False

    with open(filename, 'r') as f:
        for line in f:
            line_list = line.split()
            #print("line list", line_list)
            if line_list[0] == "1":
                line_list[0] = "Si"
            elif line_list[0] == "2":
                line_list[0] = "C"
            if line_list[0] == "3":
                line_list[0] = "Si"
            elif line_list[0] == "4":
                line_list[0] = "C"

            line = " ".join(line_list) + "\n"
            outfile.write(line)
            #print("new line", line)
    outfile.close()      
    return file_ase