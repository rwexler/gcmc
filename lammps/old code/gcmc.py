#!/root/miniconda3/bin/python -i
# preceding line should have path for Python on your machine

# mc.py
# Purpose: mimic operation of example/MC/in.mc via Python
# Syntax:  mc.py in.mc
#          in.mc = LAMMPS input script

from __future__ import print_function
import sys,random,math
import numpy as np

# set these parameters
# make sure neigh skin (in in.mc) > 2*deltamove

kb = 8.6173303e-5			# ev / k

nloop = 3000
deltaperturb = 0.5
deltamove = 0.1
T  = 500
kT = kb*T
random.seed(27848)

# parse command line

argv = sys.argv
if len(argv) != 2:
  print("Syntax: mc.py in.mc")
  sys.exit()

infile = sys.argv[1]

from lammps import lammps, LAMMPS_INT
lmp = lammps()

mu = -0.1
lmp.command("variable mu internal -0.1")

# just sets up MC problem
# run infile all at once
lmp.file(infile)

lmp.command("variable e equal pe")

# run 0 to get energy of perfect lattice
# emin = minimum energy

lmp.command("run 0")

natoms = lmp.extract_global("natoms",0)
emin = lmp.extract_compute("thermo_pe",0,0) 
lmp.command("variable emin equal $e")
fmin = emin - mu*natoms
f = emin - mu*natoms
lmp.command("variable fmin equal 'e - mu*natoms'")
lmp.command("variable f equal 'e - mu*natoms'")

# disorder the system
# estart = initial energy

x = lmp.extract_atom("x",3)
types = lmp.extract_atom("type", 0)

for i in range(natoms):
  x[i][0] += deltaperturb * (2*random.random()-1)
  x[i][1] += deltaperturb * (2*random.random()-1)
  x[i][2] += deltaperturb * (2*random.random()-1)

lmp.command("variable elast equal $e")
lmp.command("variable flast equal $f")
lmp.command("thermo_style custom step v_fmin v_flast pe")
lmp.command("run 0")
x = lmp.extract_atom("x",3)
  
estart = lmp.extract_compute("thermo_pe",0,0) 
fstart = estart - mu*natoms

# loop over Monte Carlo moves
# extract x after every run, in case reneighboring changed ptr in LAMMPS

flast = fstart
naccept = 0

act_p = [0.4, 0.3, 0.3]
#  probability of [0] MOVE, [1] ADD, [2], REMOVE

act_cdf = act_p

for i in range(len(act_cdf) - 1) :
   act_cdf[i + 1] += act_cdf[i]

"""
	For making this a GCMC simulation, use the lammps commands
	def add_atoms():
		old_natoms = natoms
		randomly choose type, x, y, z
		LOOP type, x, y, z:
			create_atoms type single x y z
		group proposed_addition id old_natoms:natoms
		
		calculate energy
		
		if REJECT:
			delete_atoms proposed_addition
	def remove_atoms():
		randomly choose rand_1, ...
		type = type of rand_1, rand_2 ..., ID'd particles
		x, y, z = positions of rand_1, rand_2, ..., ID'd particles
		group proposed_deletion id rand_1 rand_2 ...
		delete_atoms proposed_deletion
		
		calculate energy
		
		if REJECT:
			LOOP type, x, y, z:
				create_atoms type single x y z
"""
for i in range(nloop):
  natoms = lmp.extract_global("natoms",0)
  iatom = random.randrange(0,natoms)
  act_cndt = np.random.rand()
  if act_cndt < act_cdf[0]:		# MOVE
    x0 = x[iatom][0]
    y0 = x[iatom][1]
    z0 = x[iatom][2]
    x[iatom][0] += deltamove * (2*random.random()-1)
    x[iatom][1] += deltamove * (2*random.random()-1)
    x[iatom][2] += deltamove * (2*random.random()-1)
  elif act_cndt < act_cdf[1]	:	# ADD
    lmp.command("create_atoms 1 random 1 107 NULL")
    lmp.command("group proposed_addition id %i" % (natoms+1))
  else:
    type0 = types[iatom]
    x0 = x[iatom][0]
    y0 = x[iatom][1]
    z0 = x[iatom][2]
    lmp.command("group proposed_deletion id %i" % iatom)
    
  lmp.command("run 1 pre no post no")
  x = lmp.extract_atom("x",3)
  types = lmp.extract_atom("type", 0)
  e = lmp.extract_compute("thermo_pe",0,0) 

  f = e - mu*natoms			# free energy per atom
  lmp.set_variable("f", f)
  if f <= flast:
    flast = f
    lmp.command("variable flast equal $f")
    naccept += 1
  elif random.random() <= math.exp((flast-f)/kT):
    flast = f
    lmp.command("variable flast equal $f")
    naccept += 1
  else:
    if act_cndt < act_cdf[0]:		# MOVE
      x[iatom][0] = x0
      x[iatom][1] = y0
      x[iatom][2] = z0
    elif act_cndt < act_cdf[1]:		# ADD
      lmp.command("delete_atoms group proposed_addition")
    else:
      lmp.command("create_atoms %i single %f %f  %f" % (type0, x0, y0, z0))
# final energy and stats

lmp.command("variable nbuild equal nbuild")
nbuild = lmp.extract_variable("nbuild",None,0)

lmp.command("run 0")
estop = lmp.extract_compute("thermo_pe",0,0) / natoms

print("MC stats:")
print("  starting free energy =",fstart)
print("  final free energy =",fstop)
print("  minimum free energy of perfect lattice =",fmin)
print("  accepted MC moves =",naccept)
print("  neighbor list rebuilds =",nbuild)
