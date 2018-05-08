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

for i in range(100) :
	for j in range(100) :
		for k in range(100) :
			nn = bv.position_nn(xsf, (xsf.lat_vec[0]*i + xsf.lat_vec[0]*j + xsf.lat_vec[0]*k)/100)
			if nn >6 and nn < 100 : 
				print nn
				print i, j, k
	print i
