"""this module defines bond valence objects and operations"""

import numpy as np
from io import el_info
from io import xsf_info
import sys

#np.random.seed(42)

# global variable definitions

class bv(object) :
    """class for representing bv objects and operations

    parameters :
    sc_num_at : float
        number of atoms in 3x3x3 supercell
    sc_lat_vec : 3 x 3 numpy array of floats
        lattice vectors of 3x3x3 supercell
    sc_at_coord : sc_num_at x 3 numpy array of floats
        atomic coordinates for 3x3x3 supercell
    nn : numpy array of ints
        number of neighbors for coordinate(s)
    """

    def __init__(self) :
        self.sc_num_at = 0
        self.sc_lat_vec = 0
        self.sc_at_coord = 0
        self.nn = 0

    def make_sc(self, xsf) :
        """make 3x3x3 supercel"""
        range = np.array([0, -1, 1])
        self.sc_lat_vec = xsf.lat_vec * range.shape[0]
        self.sc_num_at = xsf.num_at * (range.shape[0] ** 3)
        self.sc_at_coord = np.zeros((self.sc_num_at, 3))
        start = 0
        finish = xsf.num_at
        for a_range in range :
            for b_range in range :
                for c_range in range :
                    self.sc_at_coord[start:finish, :] = xsf.at_coord + \
                        a_range * xsf.lat_vec[0, :] + \
                        b_range * xsf.lat_vec[1, :] + \
                        c_range * xsf.lat_vec[2, :]
                    start = finish
                    finish += xsf.num_at
        return self.sc_lat_vec, self.sc_at_coord

    def calc_nn(self, coord) :
        """calculate number of neighbors for coordinate(s) given"""
        rmin = 1.0
        rmax = 2.5
        self.nn = np.zeros((coord.shape[0], ))
        for i in range(coord.shape[0]) :
            nn_cnt = 0
            for j in range(self.sc_num_at) :
                bond_vec = coord[i, :] - self.sc_at_coord[j, :]
                for row in range(bond_vec.shape[0]) :
                    bond_vec[row] -= self.sc_lat_vec[row, row] * \
                        int(bond_vec[row] / self.sc_lat_vec[row, row])
                dist = np.linalg.norm(bond_vec)
                if dist > rmin and dist < rmax :
                    nn_cnt += 1
                elif 0 < dist and dist < rmin : # TQ: penalty for nonzero small value
                    nn_cnt += 100
            self.nn[i] = nn_cnt
        return self.nn
