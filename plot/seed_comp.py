import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.image as mpimg

from lmp_plotting import *
from lmp_utils import *
from lmp_analyze import *
from plot_ase import *
   
def main( filename ):
    path = "../results/"
    #filename = "muSi-5.70_muC-5.75_seed"
    filename1 = filename + "1.csv"
    filename2 = filename + "2.csv"
    filename3 = filename + "3.csv"
       
    #file_top = "images/Top "
    #file_side = "images/Side "
    imgs_top = []
    imgs_side = []
    df = {}
            
    poteng = {}
    redPot = {}
    comp = {}   
    useVMD = True
    for num in range(0, 3):
        if useVMD:
            file = path + "images/" + filename + str(num+1) + "_top.bmp"
            imgs_top.append( mpimg.imread(file) )
            file = path + "images/" + filename + str(num+1) + "_side.bmp"
            imgs_side.append( mpimg.imread(file) )
        else:
            plot_ase(path + "raw/" + filename + str(num+1) + ".xyz")
            file = path + "images/" + filename + str(num+1) + "_top_ase.png"
            imgs_top.append( mpimg.imread(file) )
            file = path + "images/" + filename + str(num+1) + "_side_ase.png"
            imgs_side.append( mpimg.imread(file) )
        file = path + "csv/" + filename + str(num+1) + ".csv"
        df[num] = pd.read_csv( file )
        
        poteng[num] = df[num]['PE'].values
        redPot[num] = df[num]['redPot'].values
        nSi = df[num]['nSi'].values 
        nC = df[num]['nC'].values 
        comp[num] = nSi/nC
        
    indices = df[0]['indx'].values
    
    # plot contents of log file
    fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(16, 12))
    for num, ax in enumerate(axs):
        # plot energies
        ax[0].plot(indices, poteng[num], 'b--', label = 'Potential')
        ax[0].set_xlabel('Step Number')
        ax[0].set_ylabel('Potential Energy (eV)', color = 'b')
        ax[0].tick_params('y', labelcolor='b')
        
        ax3 = ax[0].twinx()
        ax3.plot(indices, redPot[num], 'r--', label = 'Grand Potential Potential')
        #ax3.set_ylabel('Grand Potential Energy (eV)', color = 'r')
        ax3.tick_params('y', labelcolor='r')
            
        ax[1].plot(indices, comp[num], 'k--', linewidth=1.5, label='Si/C')
        ax[1].legend(loc='upper right')
        ax[1].set_xlabel('Step Number')
        ax[1].set_ylabel('% Composition')
        
        ax[2].imshow(imgs_top[num])
        ax[3].imshow(imgs_side[num])

    axs[0,0].set_title('Potential Energy')
    axs[0,1].set_title('Ratio of Composition nSi/nC')
    axs[0,2].set_title('Orthographic Top View')
    axs[0,3].set_title('Orthographic Side View')
    fig.tight_layout()
    if useVMD:
        plt.savefig( path + filename + "_comp.png" )
    else:
        plt.savefig( path + filename + "_comp_ase.pdf" )
        
if __name__ == "__main__":
    #files = [ "muSi-5.60_muC-5.75_seed", "muSi-5.70_muC-5.75_seed", "muSi-5.80_muC-5.85_seed",
    #           "muSi-5.90_muC-5.75_seed", "muSi-6.00_muC-5.75_seed", "muSi-5.80_muC-5.75_seed" ]
    files = [ "muSi-5.80_muC-5.65_seed" ]
    for file in files:
        for i in range(1,4): 
            tmp = file + str(i) + ".lammps"
            scrape_structs(tmp)
        main( file )