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
        
    num_elements = 2
    # right now the other method to count number of structures is too slow to actually be feasible 
    num_structs = 100001
    print("num structures:", num_structs)
    # it is less computationally expensive to make a list and append items
    # then make the data frame in one go from that list
    data = []
    
    area = 5
    # [0] Silicon, [1] Carbon
                        
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
            indx = str_ind + 1
            T = float(line_list[1])
            PE = float(line_list[2])
            redPot = float(line_list[3])
            
            nSi = int(line_list[4])
            nC = int(line_list[5])
            
            nacc_insSi = int(line_list[6])
            nacc_delSi = int(line_list[7])
            nacc_insC = int(line_list[8])
            nacc_delC = int(line_list[9])
            
            data.append( [indx, T, PE, redPot, nSi, nC, nacc_insSi, nacc_delSi, nacc_insC, nacc_delC, 0, 0, 0, 0, 0, area] )
            # Gamma0, Eslab, EbulkAgO, EbulkO, Phi
            #df.loc[i + nsurf_run, 'dir'] = dir
            
    # get contents of log file
    # for LAMMPS I could imagine I would want to plot these against time
    # Temp PotEng f_gcmcSi[4] f_gcmcSi[6] f_gcmcC[4] f_gcmcC[6] v_taccSi v_taccC 
    # essentially the potential energy and the acceptances of insert/remove/translate for Si and C respectively    
    col_names = ['indx', 'T', 'PE', 'redPot', 'nSi', 'nC', 'nacc_insSi', 'nacc_delSi', 'nacc_insC', 'nacc_delC', 
                    'GammaO', 'Eslab', 'EbulkAgO', 'EbulkO', 'Phi', 'A']
                    
    df = pd.DataFrame(data, columns = col_names, )
    df.to_csv( generate_df_filename(filename), index = False)
    return df
    
if __name__ == "__main__":
    main()