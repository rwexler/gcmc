"""this module defines bond valence objects and operations"""

import numpy as np
from io import el_info
from io import xsf_info
import sys

class bv(object) :
	"""class for representing bv objects and operations"""
	
	def __init__(self) :
		self.at_nn      = np.array([]).astype('int')   # number of neighbors for coordinate(s)
		self.r_min      = np.array([]).astype('float') # lower limit for neighbor
		self.r_max      = np.array([]).astype('float') # upper limit for neighbor
		self.lat_vec_sc = np.zeros((0,3))              # 3 x 3 x 3 lattice point from (-1,-1,-1) to (1,1,1)
		
	def copy(self) :
		"""perform explicit copy"""
		cp_self = bv()
		cp_self.at_nn      = np.array(self.at_nn)
		cp_self.r_min      = np.array(self.r_min)
		cp_self.r_max      = np.array(self.r_max)
		cp_self.lat_vec_sc = np.array(self.lat_vec_sc)
		return cp_self
	
	def init(self, xsf, el) :
		"""initialize bv"""
		self.r_min = np.array(el.r_min)
		self.r_max = np.array(el.r_max)
		self.lat_vec_sc = np.zeros((0,3))
		range = np.array([-1, 0, 1])
		for i in range :
			for j in range :
				for k in range :
					self.lat_vec_sc = np.vstack((self.lat_vec_sc, i * xsf.lat_vec[0] + j * xsf.lat_vec[1] + k * xsf.lat_vec[2]))

	def at_all_nn(self, xsf) :
		"""calculate number of neighbors for each atom"""
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

	def at_single_nn(self, xsf, at_ind, *coord) : 
		"""calculate number of neighbors for a specific atom or position"""
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
