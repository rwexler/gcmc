from ase import Atoms
from ase.visualize import view
from ase.io import read, write

from lmp_utils import convert_to_ase

def main( filename ):    
    #filename = "../results/raw/muSi-5.60_muC-5.75_seed"
    for i in range(1,4):
        tmp = filename + str(i) + ".xyz"
        plot_ase(tmp)
        #view( plot_ase(tmp) )

def plot_ase(filename):
    file_no_suf = '.'.join( filename.split('.')[:-1] )   
    #file_out = file_no_suf + ".png"
    try:
        file_ase = file_no_suf + "_ase.xyz"
        slab = read(filename, -1)
    except:
        file_ase = convert_to_ase(filename)
        slab = read(file_ase, -1)
        
    write(file_no_suf + "_top_ase.png", slab) #, rotation='10z,-80x')
    write(file_no_suf + "_side_ase.png", slab, rotation='-90x')
    return slab

if __name__ == "__main__":
    main( "../results/raw/muSi-5.60_muC-5.75_seed" )
    #main( "../results/raw/muSi-5.90_muC-5.75_seed" )
    #main( "../results/raw/muSi-6.00_muC-5.75_seed" )