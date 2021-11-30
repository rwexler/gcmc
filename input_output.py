"""this module defines input and output objects and operations"""

import numpy as np
import os
import copy

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
			f_new.write('nat = ' + str(xsf.atom_num) + ',\n')
			for line in f_old :
				if 'CELL_PARAMETERS' in line :
					f_new.write(line)
					break
				else :
					f_new.write(line)
			for row in range(3) :
				f_new.write(str(xsf.lat_vecs[row, 0]) + ' ' +
							str(xsf.lat_vecs[row, 1]) + ' ' +
							str(xsf.lat_vecs[row, 2]) + '\n')

			for line in f_old :
				if 'ATOMIC_POSITIONS' in line :
					f_new.write(line)
					break
				else :
					f_new.write(line)
			for row in range(xsf.atom_num) :
				if row in xsf.atoms_removable :
					fix_str = ' 1 1 1'
				else :
					fix_str = ' 0 0 0'
				f_new.write(el.sym[xsf.atom_type[row]] + ' ' +
							str(xsf.atom_coords[row, 0]) + ' ' +
							str(xsf.atom_coords[row, 1]) + ' ' +
							str(xsf.atom_coords[row, 2]) + fix_str + '\n')

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
	for row in range(xsf.lat_vecs.shape[0]) :
		axsf_file.write(str(xsf.lat_vecs[row, 0]) + ' ' +
						 str(xsf.lat_vecs[row, 1]) + ' ' +
						 str(xsf.lat_vecs[row, 2]) + '\n')
	return axsf_file

def upd_axsf(axsf_file, iter, xsf, el) :
	"""function for updating axsf file"""
	axsf_file.write('PRIMCOORD ' + str(iter + 1) + '\n')
	axsf_file.write(str(xsf.atom_num) + ' 1\n')
	for row in range(xsf.atom_coords.shape[0]) :
		axsf_file.write(el.sym[xsf.atom_type[row]] + ' ' +
						str(xsf.atom_coords[row, 0]) + ' ' +
						str(xsf.atom_coords[row, 1]) + ' ' +
						str(xsf.atom_coords[row, 2]) + ' ' +
						str(xsf.atom_forces[row, 0]) + ' ' +
						str(xsf.atom_forces[row, 1]) + ' ' +
						str(xsf.atom_forces[row, 2]) + '\n')
	axsf_file.flush()
