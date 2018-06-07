import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy
import sys

df = pd.read_csv('data.csv')

f, (ax1, ax2) = plt.subplots(1, 2)

ratio = df['nAg'].values / df['nO'].values
init_ratio = ratio[0]
ax1.hist(ratio, zorder = 1)
ax1.axvline(init_ratio, color = 'k', ls = '--', zorder = 2)
ax1.set_xlabel(r'$n_{\rm Ag} / n_{\rm O}$')
ax1.set_yticks([])

df['Phi'].plot(kind = 'density')
ax2.set_xlabel(r'Surface energy (J/m$^2$) at $\Delta \mu_{\rm O} = 0$ eV')
ax2.set_yticks([])
ax2.set_ylabel('')
plt.show()
