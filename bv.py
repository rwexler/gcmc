"""this module defines bond valence objects and operations"""

import numpy as np
from structure import Element_info as el_info
from structure import Structure as xsf_info
import sys

class BondValence(object) :
    """class for representing bv objects and operations

        parameters :
        sc_num_at : 
        number of atoms in 3x3x3 supercell
          sc_lat_vec : 3 x 3 numpy array of floats
        lattice vectors of 3x3x3 supercell
          sc_at_coord : sc_num_at x 3 numpy array of floats
        atomic coordinates for 3x3x3 supercell
         nn : numpy array of ints
           number of neighbors for coordinate(s)
    """
    def __init__(self, xsf = None, el = None) :
        self.at_nn      = np.array([]).astype('int')   # number of neighbors for coordinate(s)
        self.lat_vec_sc = np.zeros((0,3))              # 3 x 3 x 3 lattice point from (-1,-1,-1) to (1,1,1)
        
        if el is None:
            self.r_min      = np.array([]).astype('float') # lower limit for neighbor
            self.r_max      = np.array([]).astype('float') # upper limit for neighbor
        else:
            self.r_min = np.array(el.r_min)                 # should this have the .astype('float') to be type-safe?
            self.r_max = np.array(el.r_max)
            if xsf is not None:
                range = np.array([-1, 0, 1])
                for i in range :
                    for j in range :
                        for k in range :
                            self.lat_vec_sc = np.vstack((self.lat_vec_sc, i * xsf.lat_vec[0] + j * xsf.lat_vec[1] + k * xsf.lat_vec[2]))

    def copy(self) :
        return copy.deepcopy(self)
			
    # get number of neighbors for each atom
    def at_all_nn(self, xsf) :
        self.at_nn = np.zeros(xsf.at_num).astype('int')
        for i in range(xsf.at_num) :
            at_type = xsf.at_type[i]
            for j in range(xsf.at_num) :
                if i != j :
                    dis = np.linalg.norm(self.lat_vec_sc + (xsf.at_coord[i] - xsf.at_coord[j]), axis = 1)
                    if (dis < self.r_min[at_type]).any() :
                        self.at_nn[i] += 999
                    else :
                        self.at_nn[i] += dis[dis < self.r_max[at_type]].shape[0]
        return np.array(self.at_nn)

    # get number of neighbors for a specific atom (at certain position, or at it's original position by default)
    def at_single_nn(self, xsf, at_ind, *coord) : 
        at_nn = 0
        at_type = xsf.at_type[at_ind]
        if len(coord) == 0 : 
            coord = np.array(xsf.at_coord[at_ind])
            for i in range(xsf.at_num) :
                if i != at_ind :
                    dis = np.linalg.norm(self.lat_vec_sc + (xsf.at_coord[i] - coord), axis = 1)
                    if (dis < self.r_min[at_type]).any() :
                        at_nn += 999
                    else :
                        at_nn += dis[dis < self.r_max[at_type]].shape[0]
        return at_nn
