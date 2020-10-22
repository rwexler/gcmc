#### LOG FOR LAMMPS BRANCH OF aiGCMC PROJECT
## 10/ 15/ 20
Started up my Xenial Puppy Linux machine and downloaded a tarball for LAMMPS. I am using cmake to configure with the following flags:
-C ../cmake/presets/minimal.cmake -D BUILD_SHARED_LIBS=on -D LAMMPS_EXCEPTIONS=on \
-D PKG_PYTHON=on -D PKG_MC=on -D PKG_REPLICA=on \
-DPYTHON_LIBRARY=/root/miniconda3/lib/libpython3.7m.so \
-DPYTHON_INCLUDE_DIR=/root/miniconda3/include/python3.7m/ \
-DPYTHON_EXECUTABLE=/root/miniconda3/bin/python

cmake --install . doesn't work for me, make install does

/////////////////////////////////////////////////////////

None of that worked. It wouldn't create a shared library liblammps.so. Of course, I could create a static library liblammps.a in a 
static cmake build. Inspired by Axel, I decide to just do a traditional make build. To make sure the python install works, I replace
the Makefile.lammps with the Makefile.lammps.python3 to have it search for python3 since the Puppy also has python 2.7. 
make yes-mc
make mode=shlib serial
make install-python

I ran the python interpreter
from lammps import lammps
l = lammps()
l.command("info config")

It all worked! Next I will try to run the mc.py code and in.mc to test it all together in a real script. Unfortunately putting it on mode-shlibs
means that it does not create an executable, so I ran it again without that mode.
# to disable track pad scrolling in bash
xinput set-prop 6 "Device Enabled" 0

## 10 / 20/ 20
I found that 'fix gcmc' is by far the simplest way to implement a LAMMPS GCMC algorithm.
I'll have to figure out later how to mix that with Rob's DFT code. Maybe I'll do 'fix gcmc; run 1'
I found that lammps trajectories don't work in VMD when the number of particles change.
Use Axel's topotools command
topo readvarxyz file.xyz

## 10/ 21/ 20
I found out that by looking at the fix_gcmc.cpp, it appears that the at least the translation attemp it not a greedy Metropolis algorithm (line 856 fix_gcmc.cpp)
Since it just checks a random number against the Boltzmann weight, it does not automatically pass if the energy_after is less than energy_before.

Interestingly I ran the lammps MC for 10,000 steps which takes only around a minute, but nothing really happens to the (210) 3C-SiC lattice. It does not reconstruct as I would expect.
I figure that a large chemical potential is necessary to allow the atoms to jump to different points on the surface to reconstruct it.

It is interesting to consider that Wexler et al are interesting in surface reconstructions in which new stuff is put ontop of the bulk material. Is this still considered
"clean" surface in ultra high vacuum that are reconstructing?

The User package for fix qmmm is really more about a hybrid QM/MM in which quantum mechanical and molecular forces act on a subset of atoms. Not really what I'm looking for...

We really would just like the potential energy (maybe forces) as caculated from a self consistent field of electron waves for the surface. Rob's code makes it clear
that those parameters can be extracted from the QE outputs. I could then use FIX EXTERNAL to set additional forces on each of the atoms.
We could use this to turn off the pair wise force and have the particles evolve according to the QE forces.

"compte ID group-ID pe fix" This means that the potential energy calculation only uses the fix for the energy.

The GCMC fix full energy option means that the potential energy is called at every step, so if we wanted to only use the DFT energy for the GCMC acceptance condition
we could add the fix keyword to the potential energy compute.

We could of course make use of a hybrid QM/MM scheme in the sense that only the surface would be quantum mechanically "active" and the rest (the bulk and the gases/vacuum) will be 
driven by molecular mechanics.

## 10/ 22/ 20
I have just remembered that aiGCMC does no translations! That is why the structures look so perfect, essentially the temperature is zero, but there is a finite energy to
add or remove atoms. Also Vignesh uses VESTA to get his ball and stick image structures from an xsf.