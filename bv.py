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
		self.at_nn      = np.array([]).astype('int')
		self.lat_vec_sc = np.zeros((0,3))
		
	def copy(self) :
		cp_self = bv()
		cp_self.at_nn      = np.array(self.at_nn)
		cp_self.lat_vec_sc = np.array(self.lat_vec_sc)
		return cp_self
	
	def init(self, xsf) :
		range = np.array([0, -1, 1])
		for i in range :
			for j in range :
				for k in range :
					self.lat_vec_sc = np.vstack((self.lat_vec_sc, i * xsf.lat_vec[0] + j * xsf.lat_vec[1] + k * xsf.lat_vec[2]))
	
#	def at_all_nn(self, xsf) :
#		self.at_nn = np.zeros(xsf.at_num).astype('int')
#		for i in range(xsf.at_num) :
#			for j in range(xsf.at_num) :
