import matplotlib.pyplot as plt
import numpy as np

filename = 'coord.axsf'

# get number of structures
num_str = 0
with open(filename, 'r') as f:
    for line in f:
        if 'PRIMCOORD' in line :
            num_str += 1

# get number of atoms and
# and number of each element
# for each structure
num_at = []
comp = np.zeros((num_str, 3))
with open(filename, 'r') as f:
    for str in range(num_str) :
        for line in f:
            if 'PRIMCOORD' in line :
                break

        # get number of atoms in structure
        for line in f :
            num_at_str = int(line.split()[0])
            num_at.append(num_at_str)
            break

        # get dictionary of element counts
        el_list_str = []
        for at in range(num_at_str) :
            for line in f :
                el_list_str.append(line.split()[0])
                break
        el_dict = dict((x, el_list_str.count(x)) for x in set(el_list_str))

        # get martix with compositions
        for key, value in el_dict.items() :
            if key == 'Sr' :
                #comp[str, 0] = 1.
                comp[str, 0] = float(el_dict['Sr'])
            elif key == 'Ti' :
                comp[str, 1] = float(el_dict['Ti']) #/ float(el_dict['Sr'])
            else :
                comp[str, 2] = float(el_dict['O']) #/ float(el_dict['Sr'])

# plot composition
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(0, num_str, 1)
ax.plot(x, comp[:, 0], 'k--', linewidth=1.5, label='Sr')
ax.plot(x, comp[:, 1], 'r--', linewidth=1.5, label='Ti')
ax.plot(x, comp[:, 2], 'b--', linewidth=1.5, label='O')
ax.legend(loc='right')
ax.set_xlabel('Step Number')
ax.set_ylabel('Composition')

#plt.show()
plt.savefig('comp_plot.png')