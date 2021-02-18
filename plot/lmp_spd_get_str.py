import sys
import pandas as pd
import os
import numpy as np

# get directory containing structure of phase with user-input id
phase_id = int(sys.argv[1])
df = pd.read_csv('data.csv')
dirs = df.dir.values
str_id = df.indx.values
phase_dir = dirs[phase_id]
phase_str_id = str_id[phase_id]

# get structure
os.chdir(phase_dir)
with open('coord_new.axsf', 'r') as file :

	# get lattice vectors
	lat_vect = np.zeros((3, 3))
	for line in file :
		if 'PRIMVEC' in line :
			break
	for i in range(3) :
		for line in file :
			lat_vect[i, :] = np.array(line.split()).astype(float)
			break

	# find coordinates
	for line in file :
		if 'PRIMCOORD ' + str(phase_str_id) in line :
			break

	# get number of atoms
	for line in file :
		nat = int(line.split()[0])
		break

	# get coordinates
	el_sym = list()
	at_coord = np.zeros((nat, 3))
	for i in range(nat) :
		for line in file :
			el_sym.append(line.split()[0])
			at_coord[i, :] = np.array(line.split()[1:4]).astype(float)
			break

os.chdir('../../plot')

# make xsf file
with open(str(phase_id) + '.xsf', 'w') as file :
	file.write('CRYSTAL\nPRIMVEC\n')
	for i in range(3) :
		file.write(str(lat_vect[i, :])[1:-1] + '\n')
	file.write('PRIMCOORD\n' + str(nat) + '\n')
	for i in range(nat) :
		file.write(el_sym[i] + ' ' + str(at_coord[i, :])[1:-1] + '\n')
