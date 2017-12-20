"""this module defines input and output objects and operations"""

import numpy as np
import os
import copy

# global variable definitions
h = 6.626070040e-34 # j * s
kb = 1.38064852e-23 # j / k
amu_kg = 1.660539040e-27 # kg / amu
ry_ev = 13.605693009 # ev / ry

class xsf_info(object) :
    """class for representing an xsf file

    parameters :

    filename : str
        xsf filename
    lat_vec : 3 x 3 matrix of floats
        lattice vectors
    el_list : list of str
        list of element symbols
    num_at : int
        number of atoms
    at_coord : num_at x 3 array of floats
        atomic coordinates
    num_each_el : dict of str keys and int values
        number of each unique element
    ind_rem_at :  list of int
        list of indices of removable atoms
    c_min : float
        minimum allowed projection of atomic coordinates on c
    c_max : float
        maximum allowed projection of atomic coordinates on c
    vol : int
        volume of variable composition region
    """

    def __init__(self, filename) :
        self.filename = filename
        self.lat_vec = np.zeros((3, 3))

        with open(self.filename, 'r') as f :
            for line in f :
                if 'PRIMCOORD' in line :
                    break
            for line in f :
                self.num_at = np.asarray(line.split())[0].astype('int')
                break

        self.at_coord = np.zeros((self.num_at, 3))
        self.el_list = []
        self.num_each_el = {}
        self.ind_rem_at = []
        self.c_min = 0
        self.c_max = 0
        self.vol = 0

    def get_lat_vec(self) :
        """get lattice vectors attribute"""
        with open(self.filename, 'r') as f :
            for line in f :
                if 'PRIMVEC' in line :
                    break
            for row in range(0, 3) :
                for line in f :
                    self.lat_vec[row, :] = np.asarray(line.split()).astype('float')
                    break
        return self.lat_vec

    def get_at_coord(self) :
        """get atomic coordinates attribute"""
        with open(self.filename, 'r') as f :
            for line in f :
                if 'PRIMCOORD' in line :
                    break
            f.next()
            for at in range(self.num_at) :
                for line in f :
                    self.at_coord[at, :] = np.asarray(line.split())[1:4].astype('float')
                    break
        return self.at_coord

    def get_el_list(self) :
        """get list of elements attribute"""
        with open(self.filename, 'r') as f :
            for line in f :
                if 'PRIMCOORD' in line :
                    break
            f.next()
            for at in range(self.num_at) :
                for line in f :
                    self.el_list.append(line.split()[0])
                    break
            return self.el_list

    def get_num_each_el(self) :
        """get number of each unique element"""
        self.num_each_el = dict((x, self.el_list.count(x)) for x in set(self.el_list))
        return self.num_each_el

    def get_ind_rem_at(self) :
        """get index of removable atoms"""
        with open(self.filename, 'r') as f :
            for line in f :
                if 'PRIMCOORD' in line :
                    break
            f.next()
            for at in range(self.num_at) :
                for line in f :
                    rem = np.asarray(line.split())[4].astype('int')
                    if rem == 1 :
                        self.ind_rem_at.append(at)
                    break
        return self.ind_rem_at

    def get_c_min_max(self, buf_len) :
        """get minimum and maximum projection of atomic coordinates along c"""
        c_unit = self.lat_vec[2] / np.linalg.norm(self.lat_vec[2])
        c_proj = np.dot(self.at_coord, c_unit)
        perp = np.cross(self.lat_vec[0], self.lat_vec[1])
        perp_unit = perp / np.linalg.norm(perp)
        self.c_min = np.min(c_proj) - buf_len / np.dot(c_unit, perp_unit)
        self.c_max = np.max(c_proj) + buf_len / np.dot(c_unit, perp_unit)
        return self.c_min, self.c_max

    def get_vol(self) :
        """get volume of variable composition region"""
        c_unit = self.lat_vec[2] / np.linalg.norm(self.lat_vec[2])
        self.vol = np.dot(np.cross(self.lat_vec[0], self.lat_vec[1]), (self.c_max - self.c_min) * c_unit)
        return self.vol

    def pop_attr(self, buf_len) :
        """populate attributes"""
        self.get_lat_vec()
        self.get_at_coord()
        self.get_el_list()
        self.get_num_each_el()
        self.get_c_min_max(buf_len)
        self.get_vol()
        self.get_ind_rem_at()

    def copy(self, xsf_new) :
        """copy attributes from xsf_new to this xsf"""
        self.at_coord = copy.copy(xsf_new.at_coord)
        self.ind_rem_at = copy.copy(xsf_new.ind_rem_at)
        self.el_list = copy.copy(xsf_new.el_list)
        self.num_each_el = copy.copy(xsf_new.num_each_el)
        self.num_at = copy.copy(xsf_new.num_at)


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
        self.forces = []

    def get_final_en(self) :
        """get final energy attribute"""
        with open(self.filename, 'r') as f :
            for line in f :
                if '!' in line :
                    self.final_en = np.asarray(line.split())[4].astype('float')
                    break
        return self.final_en

    def get_forces(self) :
        """get forces attribute"""
        with open(self.filename, 'r') as f :
            for line in f :
                if '  force =' in line :
                    self.forces.append(line.split()[6:9])
        self.forces = np.array(self.forces).astype(float)
        return self.forces

def make_qe_in(filename, xsf) :
    """function for making QE input file"""
    with open(filename, 'w') as f_new :
        with open('./../templates/' + filename, 'r') as f_old :
            for line in f_old :
                if 'ibrav' in line :
                    f_new.write(line)
                    break
                else :
                    f_new.write(line)
            f_new.write('nat = ' + str(xsf.num_at) + ',\n')
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
            for row in range(xsf.num_at) :
                f_new.write(xsf.el_list[row] + ' ' +
                            str(xsf.at_coord[row, 0]) + ' ' +
                            str(xsf.at_coord[row, 1]) + ' ' +
                            str(xsf.at_coord[row, 2]) + '\n')

def init_log(filename) :
    """function for initializing log file"""
    log_file = open(filename, 'w')
    log_file.write('i | new_en | en_curr | en_low | acc | pace\n')
    return log_file

def upd_log(log_file, iter, new_en, mc_test) :
    """function for updating log file"""
    log_file.write(str('{0:4d}'.format(iter + 1)) + ' ' +
                   str('{0:.5f}'.format(new_en)) + ' ' +
                   str('{0:.5f}'.format(mc_test.g_curr)) + '  ' +
                   str('{0:.5f}'.format(mc_test.g_low)) + ' ' +
                   str('{0:2d}'.format(mc_test.acc)) + ' ' +
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

def upd_axsf(axsf_file, iter, xsf, forces) :
    """function for updating axsf file"""
    axsf_file.write('PRIMCOORD ' + str(iter + 1) + '\n')
    axsf_file.write(str(xsf.num_at) + ' 1\n')
    for row in range(xsf.at_coord.shape[0]) :
        axsf_file.write(xsf.el_list[row] + ' ' +
                        str(xsf.at_coord[row, 0]) + ' ' +
                        str(xsf.at_coord[row, 1]) + ' ' +
                        str(xsf.at_coord[row, 2]) + ' ' +
                        str(forces[row, 0]) + ' ' +
                        str(forces[row, 1]) + ' ' +
                        str(forces[row, 2]) + '\n')
    axsf_file.flush()

class el_info(object) :
    """class for representing element information file

    parameters :

    filename : str
        csv filename
    num_el : int
        number of elements
    el_sym : list of str
        list of element symbols
    at_wt : numpy array of float
        atomic weights
    ind_to_el_dict : dict of int keys, and str and float values
        maps element index to symbol, atomic weight, and thermal de broglie wavelength
    el_to_ind_dict : dict of str keys, and int and float values
        maps element symbol to index, atomic weight, and thermal de broglie wavelength
    """

    def __init__(self, filename) :
        self.filename = filename

        self.num_el = 0
        with open(self.filename, 'r') as f :
            for line in f :
                self.num_el += 1

        self.el_sym = []
        self.at_wt = np.zeros((self.num_el, ))
        self.therm_db = np.zeros((self.num_el, ))
        self.ind_to_el_dict = {}
        self.el_to_ind_dict = {}

    def get_el_sym(self) :
        """function for getting element symbols from csv file"""
        with open(self.filename, 'r') as f :
            for line in f :
                self.el_sym.append(line.split()[1])
        return self.el_sym

    def get_at_wt(self) :
        """function for getting atomic weights from csv file"""
        with open(self.filename, 'r') as f :
            for row, line in enumerate(f) :
                self.at_wt[row] = np.asarray(line.split())[2].astype('float') * amu_kg
        return self.at_wt

    def get_therm_db(self, T) :
        """function for getting thermal de broglie wavelength"""
        self.therm_db = np.sqrt(h ** 2 / (2 * np.pi * self.at_wt * kb * T)) * 1e10
        return self.therm_db

    def get_ind_to_el_dict(self) :
        """function for making dictionary of
        element index keys and
        element symbol, atomic weight, and thermal de broglie wavelenght values
        """
        for i in range(self.num_el) :
            self.ind_to_el_dict[i] = {'el_sym' : self.el_sym[i],
                                      'at_wt' : self.at_wt[i],
                                      'therm_db' : self.therm_db[i]}
        return self.ind_to_el_dict

    def get_el_to_ind_dict(self) :
        """function for making dictionary of
        element symbol keys and
        element index, atomic weight, and thermal de broglie wavelength values
        """
        for i in range(self.num_el) :
            self.el_to_ind_dict[self.el_sym[i]] = {'el_ind' : i,
                                                   'at_wt' : self.at_wt[i],
                                                   'therm_db' : self.therm_db[i]}
        return self.el_to_ind_dict

    def pop_attr(self, T_exc) :
        """populate attributes"""
        self.get_el_sym()
        self.get_at_wt()
        self.get_therm_db(T_exc)
        self.get_ind_to_el_dict()
        self.get_el_to_ind_dict()
