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
<<<<<<< HEAD
means that it does not create an executable, so I ran it again without that mode. To disable track pad scrolling in bash
xinput set-prop 6 "Device Enabled" 0

## 10/ 19/ 20
So I tried the mc.py python code with the preset in.mc, which perturbs a hexagonal lattice. This works well, so I took my old in.sic script
and got rid of the extra stuff to make it work as the input for mc.py. I need to have MANYBODY installed to use the tersoff potential.
Also this time I added
make yes-manybody
make yes-python
It seems that lammps and python work even if the python package isn't installed, so I'll see if installing it affects anything even.
Ok the installation didn't work with 'make yes-python'. Reverted it to no, and it works
I have decided that having 'mc.uvt_new_structure()' and 'mc.uvt_new_structure()' is confusing and frustrating to switch all calls of the 
method from one to another. I think it is simpler to have a child class MC_NP which overrides the method, so the calls are the same.
=======
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

## 10/ 27/ 20
I realized that really isn't detailed balance with insertions and deletions because insertions are equally placed proportional to 1/V while deletions are 1/N.
I can't really do anything about that so I just set the probability of doing a deletion to be 0.001 i.e. 0.1% so for every deletion there are 999 insertions, that should help offset 
the boon deletions had near chemical potential equilibrium.

pbc set {17.44 17.44 17.44} -all
This command combined with the gui allows one to set periodic images of the simulation.
I got a pretty simulation of a surface with sterically correct added monolayer of inhomogeneous SiC. That was with mu_si_ex  = -5 and mu_c_ex = -11.69

The following results are from over the weekend when I was messing around with tuning the chemical potential in in.bulk.gcmc.sic. I am putting these here for posterity.
-11.25 means that the minimal energy is to pull out all the atoms
-10.5 and -10 will eventually pull out all atoms
-9.75 atom number doesn't change
average excess mu_si = -10.9
average excess mu_c = -11.69
# M the average number of MC moves should be the on the order of the number of atoms in cell
# so that on an average cycle each atom moves per MC cycle
# full_energy is needed to include the total energy calculation since it is a many-body pair style

Today I finally did what Rob originally recommended, to take many npt simulations with different particle numbers and compare their potential energies.
I first started with 256 and 256 atoms and actually just used the porosity command to slowly remove atoms and then requilibrate for 1000 steps.

This meant I could analyze the potential energy for variable nSi and nC and create difference quotients for the chemical potential.
That worked fine and I got answer in the vicinity I expected. -7.9 for Si and -10 for C.

I tried Widom Insertion Method with 2048 insert attempts and I slowly made the system more porous. Interestingly the chemical potential for Si is 20 (very high since it is hard for SIlicon
to go in the interstitial spaces) but for C it is ~ 0 eV (Carbon interstices are a common defect in SiC) and as I decrease the number of atoms, 
the chemical potentials stabilize to be -8 to -6 eV for Si and -10 to -7 eV for C. I should plot that.

## 12/ 03/ 20
Back into the project. I cleaned up some of my code, simplifying my lammps scripts into just two: a GCMC one and a Widom one. Then I read the structural data from a lammps file.
I realized that I was calculating the ideal chemical potential from the TOTAL pressure i.e. 1.0 atmosphere. First of all, this is just too large. Second it should be the partial pressure of 
each of the gases. I changed myPress to be two variables pressSi and pressC, which hold the partial pressure of each species, currently set to 0.1 atm.

After simplifying my lammps scripts, I ran in.widom.sic with a surface 100 structure and then a bulk 100 structure. I saw a difference immediately! The excess chemical potential
for silicon and carbon are around -5.5 eV for the surface. For the bulk mu_Si_ex is around +20 eV and mu_C_ex is aroun -0.5 eV. This is agrees with what I recorded in October.

Also Nat and I were talking about the idea of burn in times. If I have a system start in a ground state, how likely is it that it will explore very different systems? Not likely at all, because
it will have to climb out of the funnel, especially relevant since all my systems start in a perfect lattice, which is the local ground state.

I equilibrate the 100 surface for 1000 steps at 500 K and it looks like it kind of increases the chemical potential. For Silicon it can get to -3.26 eV and Carbon's gets down to -4.4 eV.
It's variation is also higher though because the Carbon one can get into the -6.0 eV range as well.
Also all of today's experiments were with 2048 insertions, no deletions and nvt equilibration.

I ran some more experiments with higher Si or C chemical potential, still adding up to what I think is the excess surface bulk energy of SiC, -10.5 eV. The surface is surprisingly stable
but for long times I see the number of insertions and deletions of both increase, without destabilizing the surface. This is of course without nvt. 
When I did fix atom swap with higher propensity for carbon, the entire region was replaced with carbon! With similar high Silicon the same doesn't happen. Instead the carbon surfaces,
there was more carbon to start with, are removed to have a Silicon surface. Still though sometimes the Carbon monolayer would try to reform.
>>>>>>> 34a7daa4bfc6aaff2a6dcafed5797270d44eddea
