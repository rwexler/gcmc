from ase import Atoms
from ase.visualize import view
from ase.io import read, write

from lmp_utils import convert_to_ase

def main():    
    filename = "muSi-5.80_muC-5.75_seed3_lastframe.xyz"
    view( plot_ase(filename) )

def plot_ase(filename):
    file_no_suf = '.'.join( filename.split('.')[:-1] )   
    #file_out = file_no_suf + ".png"
    try:
        file_ase = file_no_suf + "_ase.xyz"
        slab = read(filename, -1)
    except:
        file_ase = convert_to_ase(filename)
        slab = read(file_ase, -1)
        
    write(file_no_suf + "_top.png", slab) #, rotation='10z,-80x')
    write(file_no_suf + "_side.png", slab, rotation='-90x')
    return slab

if __name__ == "__main__":
    main()