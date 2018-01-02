"""this module defines bond valence objects and operations"""

import numpy as np
from io import el_info
from io import xsf_info
import sys

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
		self.at_nn      = np.array([]).astype('int')   # number of neighbors for coordinate(s)
		self.r_min      = 0.0                          # lower limit for neighbor
		self.r_max      = 0.0                          # upper limit for neighbor
		self.lat_vec_sc = np.zeros((0,3))              # 3 x 3 x 3 lattice point from (-1,-1,-1) to (1,1,1)

	# enable explicit copy
	def copy(self) :
		cp_self = bv()
		cp_self.at_nn      = np.array(self.at_nn)
		cp_self.r_min      = self.r_min
		cp_self.r_max      = self.r_max
		cp_self.lat_vec_sc = np.array(self.lat_vec_sc)
		return cp_self
	
	def init(self, xsf) :
		self.r_min = 1.0
		self.r_max = 2.5
		self.lat_vec_sc = np.zeros((0,3))
		range = np.array([0, -1, 1])
		for i in range :
			for j in range :
				for k in range :
					self.lat_vec_sc = np.vstack((self.lat_vec_sc, i * xsf.lat_vec[0] + j * xsf.lat_vec[1] + k * xsf.lat_vec[2]))
	# get number of neighbors for each atom
	def at_all_nn(self, xsf) :
		self.at_nn = np.zeros(xsf.at_num).astype('int')
		for i in range(xsf.at_num) :
			for j in range(xsf.at_num) :
				if i != j :
					dis = np.linalg.norm(self.lat_vec_sc + (xsf.at_coord[i] - xsf.at_coord[j]), axis = 1)
					if (dis < self.r_min).any() :
						self.at_nn[i] += 100
					else :
						self.at_nn[i] += dis[dis < self.r_max].shape[0]
		return np.array(self.at_nn)

	# get number of neighbors for a specific atom (at certain position, or at it's original position by default)
	def at_single_nn(self, xsf, at_ind, *coord) : 
		at_nn = 0
		if len(coord) == 0 : 
			coord = np.array(xsf.at_coord[at_ind])
		for i in range(xsf.at_num) :
			if i != at_ind :
				dis = np.linalg.norm(self.lat_vec_sc + (xsf.at_coord[i] - coord), axis = 1)
				if (dis < self.r_min).any() :
					at_nn += 100
				else :
					at_nn += dis[dis < self.r_max].shape[0]
		return at_nn

	# get number of neighbors for a specific position
	def position_nn(self, xsf, coord) :
		at_nn = 0
		for i in range(xsf.at_num) :
			dis = np.linalg.norm(self.lat_vec_sc + (xsf.at_coord[i] - coord), axis = 1)
			if (dis < self.r_min).any() :
				at_nn += 100
			else :
				at_nn += dis[dis < self.r_max].shape[0]
		return at_nn


""" could try to (1) define differen r_min, r_max for different elements. (2) add tolerance to nn"""
