import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.image as mpimg

from lmp_plotting import *
from lmp_utils import *
from plot_ase import *
   
def main():
    filename = "muSi-5.80_muC-5.75_seed"
    filename1 = filename + "1.csv"
    filename2 = filename + "2.csv"
    filename3 = filename + "3.csv"
       
    file_top = "images/Top "
    file_side = "images/Side "
    imgs_top = []
    imgs_side = []
    df = {}
            
    poteng = {}
    redPot = {}
    comp = {}   
    
    for num in range(0, 3):
        if False:
            file = file_top + str(num+1) + ".png"
            imgs_top.append( mpimg.imread(file) )
            file = file_side + str(num+1) + ".png"
            imgs_side.append( mpimg.imread(file) )
        else:
            plot_ase(filename + str(num+1) + ".xyz")
            file = filename + str(num+1) + "_top.png"
            imgs_top.append( mpimg.imread(file) )
            file = filename + str(num+1) + "_side.png"
            imgs_side.append( mpimg.imread(file) )
        file = filename + str(num+1) + ".csv"
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
    plt.savefig( "seed_comp_ase.pdf" )
        
if __name__ == "__main__":
    main()