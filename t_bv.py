from bv import bv
from io import xsf_info
from io import el_info

file1 = "el_list.txt"
file2 = "structure.xsf"

el = el_info()
xsf = xsf_info()
bv = bv()

el.init(file1)
xsf.init(file2, el)
bv.init(xsf)

print bv.lat_vec_sc
