"""this module defines bond valence objects and operations"""

import numpy as np
from io import el_info
from io import xsf_info
import sys

#np.random.seed(42)

T_exc = 1
buf_len = 2.0
el = el_info('el_list.txt')
el.pop_attr(T_exc)
xsf = xsf_info('structure.xsf')
xsf.pop_attr(buf_len, el)

def supercell(xsf) :
    bound = np.array([0, -1, 1])
    sc_lat_vec = xsf.lat_vec * bound.shape[0]
    sc_num_at = xsf.num_at * (bound.shape[0] ** 3)
    sc_at_coord = np.zeros((sc_num_at, 3))
    start = 0
    finish = xsf.num_at
    for a_bound in bound :
        for b_bound in bound :
            for c_bound in bound :
                sc_at_coord[start:finish, :] = xsf.at_coord + \
                    a_bound * xsf.lat_vec[0, :] + \
                    b_bound * xsf.lat_vec[1, :] + \
                    c_bound * xsf.lat_vec[2, :]
                start = finish
                finish += xsf.num_at
    return sc_lat_vec, sc_at_coord

nneighbors = []
rcut = 2.5 # a
sc_lat_vec, sc_at_coord = supercell(xsf)
for i in range(xsf.num_at) :
    temp = 0
    for j in range(sc_at_coord.shape[0]) :
        if i != j :
            bond_vec = sc_at_coord[i, :] - sc_at_coord[j, :]
            for row in range(bond_vec.shape[0]) :
                bond_vec[row] -= sc_lat_vec[row, row] * int(bond_vec[row] / sc_lat_vec[row, row])
            dist = np.linalg.norm(bond_vec)
            if dist < rcut :
                temp += 1 
    nneighbors.append(temp)

for i in range(len(xsf.el_list)) :
    print i + 1, xsf.el_list[i], nneighbors[i]

def gen_test_pos(xsf) :
    test_pos = np.zeros((3, ))
    for row in xsf.lat_vec :
        test_pos += np.random.rand() * row
    return test_pos

def calc_nneighbors(test_pos, sc_lat_vec, sc_at_coord) :
    nn = 0
    for j in range(sc_at_coord.shape[0]) :
        bond_vec = test_pos - sc_at_coord[j, :]
        for row in range(bond_vec.shape[0]) :
            bond_vec[row] -= sc_lat_vec[row, row] * int(bond_vec[row] / sc_lat_vec[row, row])
        dist = np.linalg.norm(bond_vec)
        if dist < rcut :
            nn += 1
    return nn

niter = 1000
nneighbors2 = []
for i in range(niter) :
    test_pos = gen_test_pos(xsf)
    nneighbors2.append(calc_nneighbors(test_pos, sc_lat_vec, sc_at_coord))

print nneighbors2
print np.array(nneighbors2).min()
print np.array(nneighbors2).max()

# preferred number of neighbors for each element

print el.pref_coord
