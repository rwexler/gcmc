import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

filename = 'log.dat'

# get contents of log file
new_en = []
en_curr = []
en_low = []
acc = []
pace = []
with open(filename, 'r') as f :
    for line in f :
        if 'new_en' in line :
            break
    for line in f :
        new_en.append(float(line.split()[1]))
        en_curr.append(float(line.split()[2]))
        en_low.append(float(line.split()[3]))
        acc.append(int(line.split()[4]))
        pace.append(float(line.split()[5]))

new_en = np.asarray(new_en)
en_curr = np.asarray(en_curr)
en_low = np.asarray(en_low)
acc = np.asarray(acc)
pace = np.asarray(pace)

# plot contents of log file
fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8, 4))
x = np.arange(0, len(new_en), 1)

# plot energies
#ax0.plot(x, new_en,  'k--', label = 'New')
ax0.plot(x, en_curr, 'r--', label = 'Current')
ax0.plot(x, en_low,  'b--', label = 'Low')
ax0.legend(loc="upper right")
ax0.set_xlabel('Step Number')
en_min = np.min(en_curr)
en_max = np.max(en_curr)
ax0.set_ylabel('Free Energy (eV)')
ax0.set_ylim(en_min, en_max)

# plot pace
ax1.plot(x, acc, 'k--')
ax1.set_xlabel('Step Number')
ax1.set_ylabel('Number of Accepted Steps', color = 'k')
ax1.tick_params('y', color='k')

ax2 = ax1.twinx()
ax2.plot(x, pace, 'r--')
ax2.set_ylabel('Percent of Steps Accepted', color = 'r')
ax2.tick_params('y', color='r')

fig.tight_layout()
plt.savefig('log_plot.pdf')