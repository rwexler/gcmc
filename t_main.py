from new_io import xsf_info, el_info
from new_io import make_qe_in
import os
import sys
import numpy as np

xsf_filename = sys.argv[1]
el_filename = sys.argv[2]
T = 1

el = el_info()
xsf = xsf_info()
el.pop_attr(el_filename,T)
xsf.pop_attr(xsf_filename, el, 2.0)

print xsf.at_num
print xsf.lat_vec
print xsf.at_type
print xsf.el_each_num
print xsf.at_coord
print xsf.at_rmb
print xsf.at_swap
print xsf.c_min
print xsf.c_max
print xsf.vol

print el.num
print el.sym
print el.wt
print el.pref_coord
print el.therm_db

os.system('mkdir -p temp')
os.chdir('temp')
make_qe_in('qe.in', xsf, el)
