import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('data.csv')
dmuOs = np.linspace(-2, 0, 101)
pref_phase = np.zeros((dmuOs.shape[0], ))

for index, row in df.iterrows() :

	# calculate gammaO
	gammaO = row.loc['GammaO']
	
	# calculate phi
	phi = row.loc['Phi']
    
	# get surface area
	A = row.loc['A']

	omega = np.zeros((dmuOs.shape[0], ))
	for i, dmuO in enumerate(dmuOs) :
        
		# calculate omega in J/m^2
		omega[i] = (1. / (2. * A)) * (phi + gammaO * dmuO) * 16.0218
		# 1 eV/A^2 = 16.0218 J/m^2

	# check if preferred phase
	if index == 0 :
		omega_min = omega
	else :
		for i in range(dmuOs.shape[0]) :
			if omega[i] < omega_min[i] :
				pref_phase[i] = index
				omega_min[i] = omega[i]

	# plot surface energy line
	plt.plot(dmuOs, omega, alpha = 0.6, zorder = 1)

	# plot first surface energy line
	if index == 0 :
		plt.plot(dmuOs, omega, 'b', lw = 2, zorder = 2)

# plot minimum surface energy and print preferred phases
plt.plot(dmuOs, omega_min, 'k', lw = 2, zorder = 2)
for i in range(dmuOs.shape[0]) :
	print(dmuOs[i], pref_phase[i], omega_min[i])

# plot bulk stability region
min = -3
max = 1
xmin = np.repeat(min, 10)
xmax = np.repeat(max, 10)
y = np.linspace(0, 5, 10) 
plt.plot(xmin, y, 'k:', zorder = 3)
plt.plot(xmax, y, 'k:', zorder = 3)

# plot format
plt.xlim(-2, 0)
plt.ylim(0.8, 1.1)
plt.xlabel(r'$\Delta \mu_{\rm O}$ (eV)')
plt.ylabel(r'Surface energy (J/m$^2$)')

plt.show()
#plt.savefig('spd.pdf')
