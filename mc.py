"""this module defines monte carlo objects and operations """

import random
import copy
import numpy as np
from io import xsf_info
from io import el_info

# set random seeds
#rand_seed = 52
#random.seed(a = rand_seed)
#np.random.seed(seed = rand_seed)

# global variable definitions
ry_ev = 13.605693009
kb = 8.6173303e-5 # ev / k

class mc :
	"""class for representing monte carlo operations

	parameters :

	T_mov : int
	    temperature for moving steps
	T_exc : int
	    temperature for atom exchange steps
	etc...
	"""

	def __init__(self, T, T2, pace, xsf) :
		self.T = T                # System temperature # RBW: change to T_move
		self.T_exc = T2           # Temperature for exchange atoms # RBW : change to T_exc
		self.pace = pace          # Max displacement # RBW : change to curr_max_disp
		self.max_pace = 0.3       # Max displacement limit # RBW : change to max_disp
		self.nvt_run_cnt = 1      # Canonical ensemble running count
		self.uvt_run_cnt = 1      # Grand canonical ensemble running count
		self.check_acc = 25       # Max displacement update rate
		self.acc = 0              # Number of acceptance
		self.mv_num = 0           # Number of atoms to move
		self.mv_ind = 0           # Indices of atoms to move
		self.mv_vec = 0           # Displacement of atoms
		self.old_coord = 0        # Coordinates in last iteration
		self.new_coord = 0        # Coordinates in this iteration
		self.ad_vec = 0           # Coordinate of the added atom
		self.en_curr = 100          # Energy in this (previous) iteration
		self.g_curr = 100           # Free energy in this (previous) iteration
		self.en_low = 100           # Lowest energy
		self.g_low = 100            # Lowest free energy?
		self.coord_opt = 0        # Coordinates associated with the lowest energy
		self.xsf_opt_at_coord = copy.copy(xsf.at_coord)
		self.xsf_opt_ind_rem_at = copy.copy(xsf.ind_rem_at)
		self.xsf_opt_el_list = copy.copy(xsf.el_list)
		self.xsf_opt_num_each_el = copy.copy(xsf.num_each_el)
		self.xsf_opt_num_at = xsf.num_at
		self.uvt_rm_ind = 0       # Grand canonical ensemble, index of the atom to be removed
		self.uvt_ad_ind = 0       # Grand canonical ensemble, index of the atom to be added
		self.uvt_act = 0          # Grand canonical ensemble actions, 0: move atom, 1: add atom, -1: del atom
		self.uvt_exc_el = 0       # Grand canonical ensemble, index of the element to be exchanged

	# determine atoms to mv and step 
	def rand_mv(self, xsf) : # xsf is of xsf_info class
		self.mv_num = np.random.randint(1, xsf.num_at + 1)
		if self.mv_num == xsf.num_at :
			self.mv_ind = list(range(self.mv_num))
		else :
 			self.mv_ind = random.sample(range(xsf.num_at), self.mv_num)
		self.mv_vec = (np.random.rand(self.mv_num, 3) * 2 - 1) * self.pace

	# update new coords
	def new_coords(self, xsf) :
		self.rand_mv(xsf)
		self.old_coord = np.copy(xsf.at_coord)
		self.new_coord = np.copy(xsf.at_coord)
		for i, ind in enumerate(self.mv_ind) :
			self.new_coord[ind, :] += self.mv_vec[i, :]
	
	# update pace
	def update_pace(self) :
		"""get 
		"""
		if float(self.acc) / self.check_acc < 0.2 :
			self.pace /= 2.0
		elif float(self.acc) / self.check_acc > 0.4 :
			self.pace *= 1.618

		if self.pace > self.max_pace :
			self.pace = self.max_pace
		self.acc = 0
	
	# uvt update structure
	"""
		remove / add one atom or move atoms at each iteration
		atoms could only be added to / removed from the top of 
		the initial geometry.
	"""
	def uvt_new_structure(self,xsf,el,act_p) : # el is of el_info class
		# act_p defines probablity of taking different actions, [0]: move, [1]: swap, [2]: add, [3]: remove

		if len(xsf.ind_rem_at) == 0 : # if no atom is removable, set p_remove = 0
			act_p[3] = 0
		# normalize act_p, and make it accumalate probablity
		act_p = act_p / np.sum(act_p)
		act_p[1] += act_p[0]
		act_p[2] += act_p[1]
		act_p[3] += act_p[2]
		# generate a random number between 0 and 1
		print act_p
		exit()
		cndt = np.random.rand()
		if cndt < act_p[0] :    # move atoms
			self.uvt_act = 0
			self.uvt_exc_el = 0
			self.rand_mv(xsf)
			for i, ind in enumerate(self.mv_ind) :
				xsf.at_coord[ind, :] += self.mv_vec[i, :]
			return xsf.at_coord, xsf.ind_rem_at, xsf.el_list, xsf.num_each_el, xsf.num_at
		elif cndt < act_p[1] : # swap atoms
			return 0
		elif cndt < act_p[2] : # add one atom
			dis = 0
			trial = 0
                        self.uvt_act = 1
                        self.uvt_ad_ind = xsf.num_at                            # creat the index of atom to be added
                        self.uvt_exc_el = np.random.randint(el.num_el)          # find the element index
                        el_to_ad = el.ind_to_el_dict[self.uvt_exc_el]['el_sym'] # find element symbol
			while dis < 1.5 : # control atom distance
				self.ad_vec = np.zeros(3)
				self.ad_vec += np.random.rand() * xsf.lat_vec[0] 
				self.ad_vec += np.random.rand() * xsf.lat_vec[1]
				self.ad_vec += (np.random.rand() * (xsf.c_max - xsf.c_min) + xsf.c_min) / np.linalg.norm(xsf.lat_vec[2]) * xsf.lat_vec[2]
				dis = min([np.linalg.norm(self.ad_vec - xsf.at_coord[ind]) for ind in range(xsf.num_at)])
				trial += 1
				if trial >= 1000 :
					break
			xsf.at_coord = np.vstack((xsf.at_coord, self.ad_vec))
			xsf.ind_rem_at.append(self.uvt_ad_ind)                  # update the list of removable atoms
			xsf.el_list.append(el_to_ad)                            # add atom to the atom list
			xsf.num_each_el[el_to_ad] += 1                          # increase the number of that element
			xsf.num_at += 1                                         # increase the total number of atoms
			return xsf.at_coord, xsf.ind_rem_at, xsf.el_list, xsf.num_each_el, xsf.num_at
		else :             # remove one atom
			self.uvt_act = -1
			self.uvt_rm_ind = random.choice(xsf.ind_rem_at)                   # index of atom to be removed
			el_to_rm = xsf.el_list[self.uvt_rm_ind]                           # find the element symbol of the atom to be removed
			self.uvt_exc_el = el.el_to_ind_dict[el_to_rm]['el_ind']           # find the element index
			xsf.at_coord = np.delete(xsf.at_coord, (self.uvt_rm_ind), 0)      # remove the coordinates
			xsf.ind_rem_at.remove(self.uvt_rm_ind)                            # update the list of removable atoms
			xsf.ind_rem_at = [ind if ind < self.uvt_rm_ind else ind - 1 for ind in xsf.ind_rem_at]
			del xsf.el_list[self.uvt_rm_ind]                                  ##### remove the atom from atoms list # need to change index of atoms after this one!
			xsf.num_each_el[el_to_rm] -= 1                                    # decrease the number of that element
			xsf.num_at -= 1                                                   # decrease the total number of atoms
			return xsf.at_coord, xsf.ind_rem_at, xsf.el_list, xsf.num_each_el, xsf.num_at

	# uvt update structure for nano particles
	def uvt_new_structure_np(self,xsf,el) : # el is of el_info class
		max_dis = max([np.linalg.norm(xsf.at_coord[ind]) for ind in range(xsf.num_at)])
		if len(xsf.ind_rem_at) == 0 :                                   # avoid removing from void removable list
			cndt = np.random.rand() * 0.75
		elif max_dis > 3.5 :                                              ##### avoid np get in touch, parameter 5 need to be change to variable
			cndt = np.random.rand() * 0.75
			if cndt >= 0.5 :
				cndt += 0.25
		else :
			cndt = np.random.rand()
		if cndt < 0.50 :    # move atoms
			self.uvt_act = 0
			self.uvt_exc_el = 0
			self.rand_mv(xsf)
			for i, ind in enumerate(self.mv_ind) :
				xsf.at_coord[ind, :] += self.mv_vec[i, :]
			return xsf.at_coord, xsf.ind_rem_at, xsf.el_list, xsf.num_each_el, xsf.num_at
		elif cndt < 0.75 : # add one atom
			dis = 0
			while dis < 1.5 or dis > 2.5 : # control atom distance
				self.uvt_act = 1
				self.uvt_ad_ind = xsf.num_at                            # creat the index of atom to be added
				self.uvt_exc_el = np.random.randint(el.num_el)          # find the element index
				el_to_ad = el.ind_to_el_dict[self.uvt_exc_el]['el_sym'] # find element symbol 
				self.ad_vec = np.zeros(3)
				self.ad_vec += ( 2 * np.random.rand() - 1 ) * xsf.lat_vec[0] 
				self.ad_vec += ( 2 * np.random.rand() - 1 ) * np.random.rand() * xsf.lat_vec[1]
				self.ad_vec += ( 2 * np.random.rand() - 1 ) * np.random.rand() * xsf.lat_vec[2]
				dis = min([np.linalg.norm(self.ad_vec - xsf.at_coord[ind]) for ind in range(xsf.num_at)])

			xsf.at_coord = np.vstack((xsf.at_coord, self.ad_vec))
			xsf.ind_rem_at.append(self.uvt_ad_ind)                  # update the list of removable atoms
			xsf.el_list.append(el_to_ad)                            # add atom to the atom list
			xsf.num_each_el[el_to_ad] += 1                          # increase the number of that element
			xsf.num_at += 1                                         # increase the total number of atoms
			return xsf.at_coord, xsf.ind_rem_at, xsf.el_list, xsf.num_each_el, xsf.num_at
		else :             # remove one atom
			self.uvt_act = -1
			self.uvt_rm_ind = random.choice(xsf.ind_rem_at)                   # index of atom to be removed
			el_to_rm = xsf.el_list[self.uvt_rm_ind]                           # find the element symbol of the atom to be removed
			self.uvt_exc_el = el.el_to_ind_dict[el_to_rm]['el_ind']           # find the element index
			xsf.at_coord = np.delete(xsf.at_coord, (self.uvt_rm_ind), 0)      # remove the coordinates
			xsf.ind_rem_at.remove(self.uvt_rm_ind)                            # update the list of removable atoms
			xsf.ind_rem_at = [ind if ind < self.uvt_rm_ind else ind - 1 for ind in xsf.ind_rem_at]
			del xsf.el_list[self.uvt_rm_ind]                                  ##### remove the atom from atoms list # need to change index of atoms after this one!
			xsf.num_each_el[el_to_rm] -= 1                                    # decrease the number of that element
			xsf.num_at -= 1                                                   # decrease the total number of atoms
			return xsf.at_coord, xsf.ind_rem_at, xsf.el_list, xsf.num_each_el, xsf.num_at

	# canonical acceptance condition, return 1 if accepted, 0 otherwise
	# en in unit of eV
	def nvt_mc(self, en) :
		self.nvt_run_cnt += 1
		if en < self.en_curr :
			self.en_curr = en
			self.acc += 1
			if en < self.en_low :
				self.en_low = en
				self.coord_opt = np.copy(self.new_coord)
			if self.nvt_run_cnt % self.check_acc == 0 :
				self.update_pace()
			return 1
		elif np.random.uniform() < np.exp(-(en - self.en_curr) / (self.T * kb)) :
			self.en_curr = en
			self.acc += 1
			if self.nvt_run_cnt % self.check_acc == 0 :
				self.update_pace()
			return 1
		else :
			if self.nvt_run_cnt % self.check_acc == 0 :
				self.update_pace()
			return 0

	# Calculate free energy and probability threshhold of accepting steps
	# en in unit of eV
	def get_free_g_p(self, en, xsf, el, mu_list) : 
		free_g_new = en
		for i in range(el.num_el) :
			el_sym = el.ind_to_el_dict[i]['el_sym']
			free_g_new -= mu_list[i] * xsf.num_each_el[el_sym]
		exc_therm_db = el.ind_to_el_dict[self.uvt_exc_el]['therm_db']
		exc_el_num = xsf.num_each_el[el_sym] - (self.uvt_act - 1) / 2
		if self.uvt_act == 0 :
			exp_coef = np.exp(-(free_g_new - self.g_curr) / (self.T * kb))
		else :
			exp_coef = np.exp(-(free_g_new - self.g_curr) / (self.T_exc * kb))
		prob_acc = np.minimum(1, exp_coef * (1.33*3.14*3.5**3 / exc_therm_db**3 / exc_el_num )**self.uvt_act)  ##### changed volume to volume in NP cases
		return free_g_new, prob_acc

	# Adjust T_exc
	def update_T_const(self, T1, iter, period) :
		self.T = T1

	def update_T_linear(self, T1, iter, period) : ##### adjust to change temperature for moving rather than exchanging.
		self.T = T1 + float(1 - T1) / period * (iter % period)

	def update_T_exp(self, T1, iter, period) : ##### adjust to change temperature for moving rather than exchanging.
		self.T = T1 * (1 / float(T1))**( float(iter % period) / (period-1))

	def update_T_quadratic(self, T1, iter, period) : ##### adjust to change temperature for moving rather than exchanging.
		self.T = float(T1 - 1) / (period-1)**2 * (period - iter%period - 1)**2 + 1

	# Grand canonical acceptance condition, return 1 if accepted, 0 otherwise
	def uvt_mc(self, en, xsf, el, mu_list) :
		self.uvt_run_cnt += 1
		rand = np.random.rand()
		free_g, prob_acc = self.get_free_g_p(en, xsf, el, mu_list)
		if self.uvt_act == 0 : # move
			self.nvt_run_cnt += 1
			if rand <= prob_acc : # accept
				self.acc += 1
				self.g_curr = free_g
				if free_g < self.g_low : 
					self.g_low = free_g
					self.xsf_opt_at_coord = copy.copy(xsf.at_coord)
					self.xsf_opt_ind_rem_at = copy.copy(xsf.ind_rem_at)
					self.xsf_opt_el_list = copy.copy(xsf.el_list)
					self.xsf_opt_num_each_el = copy.copy(xsf.num_each_el)
					self.xsf_opt_num_at = xsf.num_at

				if self.nvt_run_cnt % self.check_acc == 0 :
					self.update_pace()
				return 1
			else :
				if self.nvt_run_cnt % self.check_acc == 0 :
					self.update_pace()
				return 0
		elif self.uvt_act == 1 or self.uvt_act == -1 : # exchange
			if rand <= prob_acc : # accept
				self.g_curr = free_g
				if free_g < self.g_low : 
					self.g_low = free_g
					self.xsf_opt_at_coord = copy.copy(xsf.at_coord)
					self.xsf_opt_ind_rem_at = copy.copy(xsf.ind_rem_at)
					self.xsf_opt_el_list = copy.copy(xsf.el_list)
					self.xsf_opt_num_each_el = copy.copy(xsf.num_each_el)
					self.xsf_opt_num_at = xsf.num_at
				return 1
			else :
				return 0
		else : 
			print 'Wrong action number! Undefined action!'
			exit()

	# Some thing might need to do: instead of fixing the action probability, separate acc into acc_a, acc_r, acc_m, and adjust action probability based on these numbers, (basicly acc_a and acc_r), i.e. if acc_a is large, meaning that adding atom is favored, we should increase the action probability of add atoms
