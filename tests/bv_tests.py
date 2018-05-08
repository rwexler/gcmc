import numpy as np
from io import el_info
from io import xsf_info
from bv import bv
import sys

#np.random.seed(42)

# test settings
T_exc = 1
buf_len = 2.0
el = el_info('el_list.txt')
el.pop_attr(T_exc)
xsf = xsf_info('structure.xsf')
xsf.pop_attr(buf_len, el)

# run tests
bvo = bv()

print '3x3x3 supercell test...'
bvo.make_sc(xsf)
print bvo.sc_num_at
print bvo.sc_lat_vec
print bvo.sc_at_coord

print 'nn for a test position (need double brackets for this code to work)...'
bvo.calc_nn(np.array([[0, 0, 10]]))
print bvo.nn

print 'nn for several test positions...'
bvo.calc_nn(np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [2, 2, 15]]))
print bvo.nn

print 'nn for atomic coordinates in xsf...'
bvo.calc_nn(xsf.at_coord)
print bvo.nn
