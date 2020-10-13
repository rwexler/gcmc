"""this module defines input and output objects and operations"""

import numpy as np
import os
import copy

# global variable definitions
h = 6.626070040e-34 # j * s
kb = 1.38064852e-23 # j / k
amu_kg = 1.660539040e-27 # kg / amu
ry_ev = 13.605693009 # ev / ry

class el_info(object) :
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

    # enable explicit copy
	def copy(self) :
		"""
        cp_self = el_info()
        cp_self.num        = self.num
		cp_self.sym        = copy.copy(self.sym)
		cp_self.wt         = np.array(self.wt)
		cp_self.therm_db   = np.array(self.therm_db)
		cp_self.pref_nn    = np.array(self.pref_nn)
		cp_self.r_min      = np.array(self.r_min)
		cp_self.r_max      = np.array(self.r_max)
        cp_self.p_add      = np.array(self.p_add)
        return cp_self
        """
        return copy.deepcopy(self)

	def update_therm_db(self, T) :
		"""function for updating thermal de broglie wavelengths"""
        if T is not 0:
            self.therm_db = np.sqrt(h ** 2 / (2 * np.pi * self.wt * kb * T)) * 1e10
        else:
            print("Temperature equal to zero! Thermal de Brogle wavelength calculation would cause division by zero")

    # deprecated the use of this function with _, not necesssary
    def _pop_attr(self, filename, T) :
		"""populate attributes"""
		self.init(filename)
		self.update_therm_db(T)

class xsf_info(object) :
	"""class for representing an xsf file"""
	def __init__(self, filename = None, el = None) : 
		self.lat_vec     = np.zeros((0, 3))               # lattice vectors
		self.at_coord    = np.zeros((0, 3))               # atomic coordinates
		self.at_force    = np.zeros((0, 3))               # forces on atoms
		self.at_type     = np.array([]).astype('int')     # indices of of element symbols
		self.el_each_num = np.array([]).astype('int')     # number of each unique element
		self.at_rmb      = np.array([]).astype('int')     # indices of removable atoms
		self.at_swap     = []                             # list of list of swappable atoms by element index
		self.at_num      = 0                              # number of atoms
		self.c_min       = 0                              # minimum allowed projection of atomic coordinates on c
		self.c_max       = 0                              # maximum allowed projection of atomic coordinates on c
		self.r_min       = 0                              # minimum allowed distance to origin
		self.r_max       = 0                              # maximum allowed distance to origin
		self.vol         = 0                              # volume of variable composition region

        if filename is not None:
            # get number of atoms
            self.set_num_atoms(self, filename)
            # get lattice vectors
            self.set_lat_vecs(self, filename)
		# initialize forces
		self.at_force = np.zeros((self.at_num, 3))
        if el is not None:    
            # get atom related attributes
            self.set_atom_attrs(self, el)        
        
        self.get_c_min_max(buf_len)
		self.get_r_min_max(buf_len)
		self.get_vol()
		self.get_vol_np()
        
    def copy(self) :
    '''
		cp_self = xsf_info()
		cp_self.lat_vec     = np.array(self.lat_vec)
		cp_self.at_coord    = np.array(self.at_coord)
		cp_self.at_force    = np.array(self.at_force)
		cp_self.at_type     = np.array(self.at_type)
		cp_self.el_each_num = np.array(self.el_each_num)
		cp_self.at_rmb      = np.array(self.at_rmb)
		cp_self.at_swap     = copy.deepcopy(self.at_swap)
		cp_self.at_num      = self.at_num
		cp_self.c_min       = self.c_min
		cp_self.c_max       = self.c_max
		cp_self.r_min       = self.r_min
		cp_self.r_max       = self.r_max
		cp_self.vol         = self.vol
		return cp_self
    '''
        return copy.deepcopy(self)

    def set_num_atoms(self, filename):
        # get number of atoms
		with open(filename, 'r') as f :
			for line in f :
				if 'PRIMCOORD' in line :
					break
			for line in f :
				self.at_num = np.array(line.split())[0].astype('int')
				break
                
    def set_lat_vecs(self, filename):
        # get lattice vectors
		self.lat_vec = np.zeros((0, 3))
		with open(filename, 'r') as f :
			for line in f :
				if 'PRIMVEC' in line :
					break
			for row in range(3) :
				for line in f :
					self.lat_vec = np.vstack((self.lat_vec, np.array([line.split()[0:3]]).astype('float')))
					break
                    
    def set_atom_attrs(self, el):
        # get atom related attributes
		self.el_each_num = np.zeros(el.num).astype('int')
		for t1 in range(el.num) :
			self.at_swap.append([])
		with open(filename, 'r') as f :
			for line in f :
				if 'PRIMCOORD' in line :
					break
			f.next()
			for at in range(self.at_num) :
				for line in f :
					# indices of element symbols of atoms
					self.at_type = np.append(self.at_type, el.sym.index(line.split()[0]))
					# number of each element
					self.el_each_num[el.sym.index(line.split()[0])] += 1
					# coordiantes of atoms
					self.at_coord = np.vstack((self.at_coord, np.array([line.split()[1:4]]).astype('float')))
					# removable atoms
					if line.split()[4] != '0' :
						self.at_rmb = np.append(self.at_rmb, at)
					# swappable atoms
					if line.split()[4] != '0' :
						self.at_swap[el.sym.index(line.split()[0])].append(at)
					break
                    
	def get_c_min_max(self, buf_len) :
		"""get minimum and maximum projection of atomic coordinates along c"""
		c_unit = self.lat_vec[2] / np.linalg.norm(self.lat_vec[2])
		c_proj = np.dot(self.at_coord, c_unit)
		perp = np.cross(self.lat_vec[0], self.lat_vec[1])
		perp_unit = perp / np.linalg.norm(perp)
		#self.c_min = np.min(c_proj) - buf_len / np.dot(c_unit, perp_unit) # exchange on both top and bottom surfaces
		#self.c_min = np.max(c_proj) - buf_len / np.dot(c_unit, perp_unit) # exchance only on top surface + a little bulk
		self.c_min = np.max(c_proj) # exchance only on top surface
		self.c_max = np.max(c_proj) + buf_len / np.dot(c_unit, perp_unit)
		return self.c_min, self.c_max
	
	def get_r_min_max(self, buf_len) :
		"""get minimum and maximum allowed r distance for NP growth"""
		r_surf = np.max(np.linalg.norm(self.at_coord, axis=1))
		self.r_min = 0.0
		self.r_max = r_surf + buf_len

	def get_vol(self) :
		"""get volume of variable composition region"""
		c_unit = self.lat_vec[2] / np.linalg.norm(self.lat_vec[2])
		self.vol = np.dot(np.cross(self.lat_vec[0], self.lat_vec[1]), (self.c_max - self.c_min) * c_unit)
		return self.vol

	def get_vol_np(self) :
		"""get volume of variable composition region for nanoparticle"""
		self.vol = 4.0 / 3 * np.pi * (self.r_max**3 - self.r_min**3)
		return self.vol
		
    # deprecated the use of this function with _, not necesssary        
	def _pop_attr(self, filename, el, buf_len) :
		"""populate attributes"""
		self.init(filename, el)
		self.get_c_min_max(buf_len)
		self.get_r_min_max(buf_len)
		self.get_vol()
		self.get_vol_np()

class qe_out_info(object) :
	"""class for representing a QE output file

	parameters :

	filename : str
		QE output filename
	final_en : float
		final energy in Ry
	forces : num_at x 3 array of floats
		forces
	"""
	
	def __init__(self, filename) :
		self.filename = filename
		self.final_en = 0
		self.forces = np.zeros((0, 3))
		self.coord = np.zeros((0, 3))

	def get_final_en(self) :
		"""get final energy attribute"""
		self.final_en = float(os.popen("grep ! " + self.filename + " | tail -1 | cut -d '=' -f 2 | cut -d 'R' -f 1").read()) * ry_ev
		return self.final_en

	def get_forces(self, at_num) :
		"""get forces attribute"""
		self.forces = np.array(os.popen("grep '  force = ' " + self.filename + " | tail -" + str(at_num) + " | cut -d '=' -f 2").read().split()).astype(float).reshape((at_num, 3))
		return np.array(self.forces)

	def get_coord(self, at_num) :
		""" get atomic coordinates attribute """

		nsteps = int(os.popen("grep ! " + self.filename + " | wc -l").read())
		# when relaxation proceeds more than one step
		if nsteps > 1 :
			command = "grep 'ATOMIC_POSITIONS' -A " + str(at_num) + " " + self.filename  + " | tail -n " + str(2*at_num+2) + " | head -n " + str(at_num) + " | sed 's/^..//g' | cut -b -50"
		# when relaxation only takes one step; possible when energy/force thresholds are large
		else :
			command = "grep 'ATOMIC_POSITIONS' -A " + str(at_num) + " " + self.filename + " | grep -v 'ATOMIC_POSITIONS' | sed 's/^..//g' | cut -b -50"
		self.coord = np.array(os.popen(command).read().split()).reshape((at_num, 3)).astype(float)
		return np.array(self.coord)

def make_qe_in(filename, xsf, el) :
	"""function for making QE input file"""
	with open(filename, 'w') as f_new :
		with open('./../templates/' + filename, 'r') as f_old :
			for line in f_old :
				if 'ibrav' in line :
					f_new.write(line)
					break
				else :
					f_new.write(line)
			f_new.write('nat = ' + str(xsf.at_num) + ',\n')
			for line in f_old :
				if 'CELL_PARAMETERS' in line :
					f_new.write(line)
					break
				else :
					f_new.write(line)
			for row in range(3) :
				f_new.write(str(xsf.lat_vec[row, 0]) + ' ' +
							str(xsf.lat_vec[row, 1]) + ' ' +
							str(xsf.lat_vec[row, 2]) + '\n')

			for line in f_old :
				if 'ATOMIC_POSITIONS' in line :
					f_new.write(line)
					break
				else :
					f_new.write(line)
			for row in range(xsf.at_num) :
				if row in xsf.at_rmb :
					fix_str = ' 1 1 1'
				else :
					fix_str = ' 0 0 0'
				f_new.write(el.sym[xsf.at_type[row]] + ' ' +
							str(xsf.at_coord[row, 0]) + ' ' +
							str(xsf.at_coord[row, 1]) + ' ' +
							str(xsf.at_coord[row, 2]) + fix_str + '\n')

def init_log(filename) :
	"""function for initializing log file"""
	log_file = open(filename, 'w')
	log_file.write('i | new_en | en_curr | en_low | acc | pace\n')
	return log_file

def upd_log(log_file, iter, new_en, mc_test) :
	"""function for updating log file"""
	log_file.write(str('{0:4d}'.format(iter + 1)) + ' ' +
				   str('{0:.5f}'.format(new_en)) + ' ' +
				   str('{0:.5f}'.format(mc_test.curr_g)) + '  ' +
				   str('{0:.5f}'.format(mc_test.opt_g)) + ' ' +
				   str('{0:2d}'.format(mc_test.nvt_acc)) + ' ' +
				   str('{0:.2f}'.format(mc_test.pace)) + '\n')
	log_file.flush()

def init_axsf(filename, niter, xsf) :
	"""function for initializing axsf file"""
	axsf_file = open(filename, 'w')
	axsf_file.write('ANIMSTEPS ' + str(niter) + '\n')
	axsf_file.write('CRYSTAL\n')
	axsf_file.write('PRIMVEC\n')
	for row in range(xsf.lat_vec.shape[0]) :
		axsf_file.write(str(xsf.lat_vec[row, 0]) + ' ' +
						 str(xsf.lat_vec[row, 1]) + ' ' +
						 str(xsf.lat_vec[row, 2]) + '\n')
	return axsf_file

def upd_axsf(axsf_file, iter, xsf, el) :
	"""function for updating axsf file"""
	axsf_file.write('PRIMCOORD ' + str(iter + 1) + '\n')
	axsf_file.write(str(xsf.at_num) + ' 1\n')
	for row in range(xsf.at_coord.shape[0]) :
		axsf_file.write(el.sym[xsf.at_type[row]] + ' ' +
						str(xsf.at_coord[row, 0]) + ' ' +
						str(xsf.at_coord[row, 1]) + ' ' +
						str(xsf.at_coord[row, 2]) + ' ' +
						str(xsf.at_force[row, 0]) + ' ' +
						str(xsf.at_force[row, 1]) + ' ' +
						str(xsf.at_force[row, 2]) + '\n')
	axsf_file.flush()
