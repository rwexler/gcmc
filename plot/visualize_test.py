from ase import Atoms
from ase.visualize import view
from ase.io import read, write

from lmp_utils import convert_to_ase

def main():
    filename = "muSi-5.80_muC-5.75_seed3_lastframe.xyz"
    file_no_suf = '.'.join( filename.split('.')[:-1] )   
                
    try:
        filename = file_no_suf + "_ase.xyz"
        slab = read(filename, -1)
    except:
        filename = convert_to_ase(filename)
        slab = read(filename, -1)
    view(slab)
    write('slab.png', slab, rotation='10z,-80x')

if __name__ == "__main__":
    main()