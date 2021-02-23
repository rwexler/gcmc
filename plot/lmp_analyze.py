import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lmp_plotting import *
from lmp_utils import *
   
def main():
    parser = argparse.ArgumentParser(
                description='Plot Composition or Statistics of LAMMPS Structure')
    parser.add_argument('--file', '-f', type=str, default='log.lammps', help='filename to scrape statistics from (default: log.lammps)')
    parser.add_argument('--scrape', dest='scrape', action='store_const', 
                const=True, default=False,
                help='scrape stats from log file instead of loading from CSV file (default: False)')
    parser.add_argument('--comp', dest='plot_comp', action='store_const', 
                const=True, default=False,
                help='plot the compositional evolution of the system (default: False)')
    parser.add_argument('--stats', dest='plot_stats', action='store_const', 
                const=True, default=False,
                help='plot the statistics of the system (default: False)')
    parser.add_argument('--acc', dest='plot_accs', action='store_const', 
                const=True, default=False,
                help='plot the acceptance numbers and percentages of the system (default: False)')
    args = parser.parse_args()
    params = vars(args)

    filename = params['file']
    if params['scrape']:
        df = scrape_structs( filename )
    else:
        df = pd.read_csv( generate_df_filename( filename ) )
    if params['plot_comp']:
        plot_comp(df, filename)
    if params['plot_stats']:
        plot_stats(df, filename)
    if params['plot_accs']:
        plot_accs(df, filename)
    
def scrape_structs(filename):
    ''' 
    Extract information on each structure from a log file into a panda DataFrame
    Save the DataFrame as a csv
    '''
    # get number of structures
    # find this line
    # Step Temp PotEng v_nSi v_nC f_gcmcSi[4] f_gcmcSi[6] f_gcmcC[4] f_gcmcC[6] v_taccSi v_taccC 
    # then count the number of lines between that and either the end of the file or a line that looks like this
    # Loop time of 1031.05 on 1 procs for 1000 steps with 505 atoms
    # the area of simulation is constant

    # WHAT WAS THE POINT OF ALL OF THESE DIAGNOSTICS? WHAT IS NEW/CURRENT/LOW ENERGY?
    new_en = []
    en_curr = []
    en_low = []
    acc = []            # Number of Accepted Steps
    pace = []           # Percent of Steps Accepted
    
    # get contents of log file
    # for LAMMPS I could imagine I would want to plot these against time
    # Temp PotEng f_gcmcSi[4] f_gcmcSi[6] f_gcmcC[4] f_gcmcC[6] v_taccSi v_taccC 
    # essentially the potential energy and the acceptances of insert/remove/translate for Si and C respectively    
    col_names = ['indx', 'T', 'PE', 'rU', 'nSi', 'nC', 'nacc_insSi', 'nacc_delSi', 'nacc_insC', 'nacc_delC', 
                    'GammaO', 'Eslab', 'EbulkAgO', 'EbulkO', 'Phi', 'A']
    num_elements = 2
    # right now the other method to count number of structures is too slow to actually be feasible 
    num_structs = 10001
    print("num structures:", num_structs)
    df = pd.DataFrame(columns = col_names, index = np.arange(num_structs))

    area = 5
    # [0] Silicon, [1] Carbon
                        
    # get number of atoms and
    # and number of each element
    # for each structure
    line_list = None
    # don't break on 'Step' because a run 0 would have a dump header with 'Step'
    # only the real run (run 1000000) would have 'PotEng' in dump header
    with open(filename, 'r') as f:
        for line in f:
            if 'PotEng' in line:
                break
        for str_ind in range(num_structs) :
            line_list = f.readline().split()
            #print("line_list:", line_list)
            #num_atoms.append( natoms_si + natoms_c )
            
            df.loc[str_ind, 'indx'] = str_ind + 1
            df.loc[str_ind, 'T'] = float(line_list[1])
            df.loc[str_ind, 'PE'] = float(line_list[2])
            df.loc[str_ind, 'redPot'] = float(line_list[3])
            
            df.loc[str_ind, 'nSi'] = int(line_list[4])
            df.loc[str_ind, 'nC'] = int(line_list[5])
            
            df.loc[str_ind, 'nacc_insSi'] = int(line_list[6])
            df.loc[str_ind, 'nacc_delSi'] = int(line_list[7])
            df.loc[str_ind, 'nacc_insC'] = int(line_list[8])
            df.loc[str_ind, 'nacc_delC'] = int(line_list[9])
            
            #df.loc[str_ind, 'GammaO'] = 0
            #df.loc[str_ind, 'Eslab'] = 0
            #df.loc[str_ind, 'EbulkAgO'] = 0
            #df.loc[str_ind, 'EbulkO'] = 0
            #df.loc[str_ind, 'Phi'] = 0
            #df.loc[str_ind, 'A'] = area
            #df.loc[i + nsurf_run, 'dir'] = dir
    df.to_csv( generate_df_filename(filename), index = False)
    return df
    
if __name__ == "__main__":
    main()