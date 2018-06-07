import pandas as pd
import glob
import os
import numpy as np

# parameters
EbulkAgO = -2407.209
EbulkO = -428.065
x = 2
y = 1

col_names = ['A', 'nAg', 'nO', 'GammaO', 'Eslab', 'EbulkAgO', 'EbulkO', 'Phi', 'dir', 'indx']
df = pd.DataFrame(columns = col_names)

nsurf_run = 0
dirs = glob.glob('../B_*/temp/')
for dir in dirs :
	os.chdir(dir)
	
	# get Gslab
	Gslab = list()
	with open('log.dat', 'r') as log :
		next(log)
		for line in log :
			Gslab.append(float(line.split()[1]))
	Gslab = np.array(Gslab)
	
	# get number of surfaces
	nsurf = len(Gslab)

	# get A, nAg, and nO
	nAg = list()
	nO = list()

	with open('coord_new.axsf', 'r') as coord :

		# get A
		bas_vec = np.zeros((2, 2))
		for line in coord :
			if 'PRIMVEC' in line :
				break
		for row in range(2) :
			for line in coord :
				bas_vec[row, :] = np.array(line.split()[0:2]).astype(float)
				break
		A = np.linalg.det(bas_vec)
		A = np.repeat(A, nsurf)

		# get nAg and nO
		for i in range(nsurf) :
			for line in coord :
				if 'PRIMCOORD' in line :
					break
			for line in coord :
				nat = int(line.split()[0])
				break
			elm = list()
			for j in range(nat) :
				for line in coord :
					elm.append(line.split()[0])
					break
			nAg.append(elm.count('Ag'))
			nO.append(elm.count('O'))
	nAg = np.array(nAg)
	nO = np.array(nO)

	# get muAg and muO
	with open('../main.py', 'r') as main :
		for line in main :
			if 'mu_list' in line :
				muAg = float(line.split()[2][1:-1])
				muO = float(line.split()[3][:-1])
				break

	# calculate gammaO
	gammaO = nAg * y / x - nO

	# calculate internal energy change
	Eslab = list()	
	for i in range(nsurf) :
		Eslab.append(Gslab[i] + nAg[i] * muAg + nO[i] * muO)
	Eslab = np.array(Eslab)

	# calculate phi
	phi = list()
	phi = Eslab - nAg * EbulkAgO / x + gammaO * EbulkO

	# make df
	for i in range(nsurf) :
		df.loc[i + nsurf_run, 'A'] = A[i]
		df.loc[i + nsurf_run, 'nAg'] = nAg[i]
		df.loc[i + nsurf_run, 'nO'] = nO[i]
		df.loc[i + nsurf_run, 'GammaO'] = gammaO[i]
		df.loc[i + nsurf_run, 'Eslab'] = Eslab[i]
		df.loc[i + nsurf_run, 'EbulkAgO'] = EbulkAgO
		df.loc[i + nsurf_run, 'EbulkO'] = EbulkO
		df.loc[i + nsurf_run, 'Phi'] = phi[i]
		df.loc[i + nsurf_run, 'dir'] = dir
		df.loc[i + nsurf_run, 'indx'] = i + 1
	nsurf_run += nsurf

	print(dir)

	os.chdir('../../plot/')

df.to_csv('data.csv', index = False)
