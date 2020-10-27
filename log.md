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

I've been experimenting with different values of the chemical potential. I was surprised to see that for anything < -6 ev for mu_c and mu_si
the surface "sublimates" or "evaporates". For > -3 eV, many more atoms are added, so I think that I want the chemical potentials
to be somewhere in between those values.

Nat and I were discussing what we are trying to do in this project: to find surface reconstructions. It is interesting to consider that when a surface
is heated up, a "eutectic" or meltly surface forms before the bulk melting point, because the free energy landscape for surface atoms will become
non-convex at lower temperatures. Could we look at the dynamics of raising the temperature to this point and cooling to see the true surface
free energy minima?

I want to see what happens when I run "minimize" for a bulk SiC lattice and amorphous bulk SiC. For the first I should expect that the tersoff potential is
at a critical point already with the correct lattice constants. For amorphous bulk SiC, will a minimize schema, such as conjugate gradient or 
dampened dynamics allow the bulk to crystallize into the correct lattice structure?

Then I want to try similar experiments with a surface of crystal and amorphous SiC. What happens when I do dynamics for elevated temperatures for these systems?

I want to talk to Rappe about how the long scale thermodynamics of a molecular dynamic system and an MC system differ. For long scales or averages, I think they will be the same.

## 10/ 23/ 20
Rob and I talked today about the way that I'm getting the wrong chemical potential, because the surface is condensing atoms rapidly. He recommends that I
read Frenkel Smit "Understanding Molecular Simulations" since they have a few chapters on Monte Carlo simulations.
He explained that elements in bulk, like silver, have a non-zero standard gibbs free energy of formation because it is elemental ideal gases, which don't interact, set the formation to zero.
It takes energy to form the bulk lattice, so that of course has an associated chemical potential that's nonzero at STP. 
I could also think about it as the vacancy chemical potential, there is an energy associated with the formation of a vacancy which is of course non zero. 
The Universitaet Kiel Professor has an interesting analogy to help illuminate the meaning of the chemical potential.
Chemical potential is the derivative of a potential (free energies) with respect to a generalized coordinate (number) so it is a generalized force.
Das macht Sinn. The equilibrium condition for chemical potentials is NOT that they are minimized ( like with the real thermodynamic potentials), instead it is
A + B -> C + D has the equilibrium condition: mu_A + mu_B = mu_C + mu_D. It the like the equilibrium of physical forces!

I should read Mathias Schaeffler's paper on Replica Exchange, he cites Rob's paper, and attempts to use ab initio dynamics for the reconstruction of surfaces, like what we are doing.

Rob recommended that a try an insertion / deletion method for an npt equilibrated lattice of SiC. That is to say measure the bulk energy per atom and see how it
changes for an identical system with one less or one more Carbon or Silicon atom. That sets up a difference quotient for a calculation of the chemical potential.

Once I get the right chemical potential, I should be able to do a lattice simulation with a high silicon chemical potential ( low carbon mu for equilibrium) and see domains of 
pure silicon form out of the SiC lattice.
## 10/26/20
I was researching this insertion method and I found it is called the Widom Insertion Method. 
The Hemholtz Free Energy is proportional to the logarithm of the partition function, A = -kT ln Q(N,V,T). Q(N,V,T) = V^N/(Lambda^(3N)*N!) * Integral{ds^N * exp[-U(s^N)/kT]},
s^N are scaled parameters. A = -kT ln V^N/(Lambda^(3N)*N!) - kT ln Integral{ds^N * exp[-U(s^N)/kT]}, so
mu = mu_ideal + mu_excess, since the first term is the ideal gas chemical potential and the excess arises from interactions.

Over the weekend I ran experiments in GCMC just to see what the "threshold" chemical potentials were for a bulk SiC lattice. That is to say, I did a binary search of mu from [-20, -5] and moved
through the interval looking for the point from which nothing happens to the lattice (because the chemical potential is too small to overcome the vacancy energy with a mild temperature of 500 K)
and when the lattice quickly sublimates (recall that is should be easier as more atoms are removed creating a cascade to zero atoms in the simulation)
Using this search I found that the critical point is between -10 eV and -9 eV. Now I realized this doesn't mean that the chemical potential of equilibrium are in that range, because even at
-10 eV it is only sometimes probable an atom is removed, but over the course of a 1000 timesteps, enough deletions occur to create a cascade. This is simply setting a 
metric for the range at which the chemical potential and vacancy formation energy on are the same order to occur with some non negible probability.

I found out that someone in July created a fix widom for LAMMPS! I tried it out and unfortnately it gave me a chemical potential of +20 eV for Si and ~ +3 eV for C. I realized it is 
because the energy of insertion is very high in a filled lattice! I needed to change the source code to add a number of attempts of deletion. I was able to stitch together
the fix widom and fix gcmc code (since fix widom is based on gcmc) and create a FixWidom::attempt_atom_deletion_full() for the full_energy calculation. I recomplied the 
code and it seems to work! I get an average chemical potential estimator of -11.15 eV for Si and -11.7 eV for C which is roughly what I expected! 
This is the excess chemical potential, while LAMMPS requires that mu = mu_ideal + mu_excess