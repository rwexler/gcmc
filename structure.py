import numpy as np
from lammps import lammps

# global variable definitions
h = 6.626070040e-34 # j * s
kb = 1.38064852e-23 # j / k
amu_kg = 1.660539040e-27 # kg / amu
ry_ev = 13.605693009 # ev / ry

class Element_info :
    """class for representing element information file"""
    def __init__(self, filename = None, T = 0) :
        self.num        = 0                              # number of elements
        self.sym        = []                             # element symbols
        self.wt         = np.array([]).astype('float')   # atomic weights
        self.therm_db   = np.array([]).astype('float')   # thermal de broglie wavelengths
        self.pref_nn    = np.array([]).astype('int')     # preferred coordination numbers
        self.r_min      = np.array([]).astype('float')   # minimum neighbor distance
        self.r_max      = np.array([]).astype('float')   # maximum neighbor distance
        self.p_add      = np.array([]).astype('float')   # probablity of choosing elements to add
        
        if filename is not None:
            # get number of elements
            with open(filename, 'r') as f :
                next(f)
                for line in f :
                    self.num += 1
            # get element info
            with open(filename, 'r') as f :
                next(f)
                for line in f :
                    # element symbol
                    self.sym.append(line.split()[1])
                    # atomic weight
                    self.wt = np.append(self.wt, np.array(line.split()[2]).astype('float') * amu_kg)
                    # preferred coordination number
                    self.pref_nn = np.append(self.pref_nn, np.array(line.split()[3]).astype('int'))
                    # minimum neighbor distance
                    self.r_min = np.append(self.r_min, np.array(line.split()[4]).astype('float'))
                    # maximum neighbor distance
                    self.r_max = np.append(self.r_max, np.array(line.split()[5]).astype('float'))
                    # probablity of choosing elements to add
            self.p_add = np.append(self.p_add, np.array(line.split()[6]).astype('float'))
            self.p_add = self.p_add / np.sum(self.p_add)
        self.update_therm_db(T)

    def copy(self) :
        return copy.deepcopy(self)

    def update_therm_db(self, T) :
        """function for updating thermal de broglie wavelengths"""
        if T is not 0:
            self.therm_db = np.sqrt(h ** 2 / (2 * np.pi * self.wt * kb * T)) * 1e10
        else:
            print("Temperature equal to zero! Thermal de Brogle wavelength calculation would cause division by zero")
    
class Structure:
    """ 
    class for representing structure of material, i.e. stores atom coordinates 
    class for representing a structure from xsf file
    """
    def __init__(self, filename = None, el = None, buf_len = 0, natoms = 0) :       
        self.atom_num      = natoms                             # number of atoms
        self.atom_coords    = np.zeros((self.atom_num, 3))               # atomic coordinates        
        self.atom_type     = np.array([]).astype('int')     # indices of of element symbols
        self.num_each_element = np.array([]).astype('int')     # number of each unique element
        self.atoms_removable      = np.array([]).astype('int')     # indices of removable atoms
        self.atoms_swap     = []                             # list of list of swappable atoms by element index
        
        self.lat_vecs     = np.zeros((0, 3))               # lattice vectors
        self.c_min       = 0                              # minimum allowed projection of atomic coordinates on c
        self.c_max       = 0                              # maximum allowed projection of atomic coordinates on c
        self.r_min       = 0                              # minimum allowed distance to origin
        self.r_max       = 0                              # maximum allowed distance to origin
        self.vol         = 0                              # volume of variable composition region
        if filename is not None:
            # get number of atoms
            self.set_num_atoms(self, filename)
            # get lattice vectors
            self.set_lat_vecss(self, filename)
            if el is not None:    
                # get atom related attributes
                self.set_atom_attrs(self, filename, el)               
            self.get_c_min_max(buf_len)
            self.get_r_min_max(buf_len)
            self.get_vol()
            #self.get_vol_np()
        
        # initialize forces on atoms
        self.atom_forces = np.zeros((self.atom_num, 3))
        
    def copy(self) :
        return copy.deepcopy(self)
        
    def set_num_atoms(self, filename):
        # get number of atoms
        with open(filename, 'r') as f :
            for line in f :
                if 'PRIMCOORD' in line :
                    break
            for line in f :
                self.atom_num = np.array(line.split())[0].astype('int')
                break
                
    def set_lat_vecs(self, filename):
        # get lattice vectors
        self.lat_vecs = np.zeros((0, 3))
        with open(filename, 'r') as f :
            for line in f :
                if 'PRIMVEC' in line :
                    break
            for row in range(3) :
                for line in f :
                    self.lat_vecs = np.vstack((self.lat_vecs, np.array([line.split()[0:3]]).astype('float')))
                    break
                    
    def set_atom_attrs(self, filename, el):
        # get atom related attributes
        self.num_each_element = np.zeros(el.num).astype('int')
        for t1 in range(el.num) :
            self.atoms_swap.append([])
        with open(filename, 'r') as f :
            for line in f :
                if 'PRIMCOORD' in line :
                    break
            f.next()
            for at in range(self.atom_num) :
                for line in f :
                    # indices of element symbols of atoms
                    self.atom_type = np.append(self.atom_type, el.sym.index(line.split()[0]))
                    # number of each element
                    self.num_each_element[el.sym.index(line.split()[0])] += 1
                    # coordiantes of atoms
                    self.atom_coords = np.vstack((self.atom_coords, np.array([line.split()[1:4]]).astype('float')))
                    # removable atoms
                    if line.split()[4] != '0' :
                        self.atoms_rm = np.append(self.atoms_rm, at)
                    # swappable atoms
                    if line.split()[4] != '0' :
                        self.atoms_swap[el.sym.index(line.split()[0])].append(at)
                    break
                    
    def get_c_min_max(self, buf_len) :
        """get minimum and maximum projection of atomic coordinates along c"""
        c_unit = self.lat_vecs[2] / np.linalg.norm(self.lat_vecs[2])
        c_proj = np.dot(self.atom_coords, c_unit)
        perp = np.cross(self.lat_vecs[0], self.lat_vecs[1])
        perp_unit = perp / np.linalg.norm(perp)
        #self.c_min = np.min(c_proj) - buf_len / np.dot(c_unit, perp_unit) # exchange on both top and bottom surfaces
        #self.c_min = np.max(c_proj) - buf_len / np.dot(c_unit, perp_unit) # exchance only on top surface + a little bulk
        self.c_min = np.max(c_proj) # exchance only on top surface
        self.c_max = np.max(c_proj) + buf_len / np.dot(c_unit, perp_unit)
        return self.c_min, self.c_max
    
    def get_r_min_max(self, buf_len) :
        """get minimum and maximum allowed r distance for NP growth"""
        r_surf = np.max(np.linalg.norm(self.atom_coords, axis=1))
        self.r_min = 0.0
        self.r_max = r_surf + buf_len

    def get_vol(self) :
        """get volume of variable composition region"""
        c_unit = self.lat_vecs[2] / np.linalg.norm(self.lat_vecs[2])
        self.vol = np.dot(np.cross(self.lat_vecs[0], self.lat_vecs[1]), (self.c_max - self.c_min) * c_unit)
        return self.vol

    def get_vol_np(self) :
        """get volume of variable composition region for nanoparticle"""
        self.vol = 4.0 / 3 * np.pi * (self.r_max**3 - self.r_min**3)
        return self.vol
