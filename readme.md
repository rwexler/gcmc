ideas

1. force-biased displacements
2. hard sphere minimum distance rule based on atomic/ionic/metallic radius
3. minimum coordination rule for nps


things to do for automated surface reconstruction

1. use to find ti phases
2. use to study nanostructure growth
3. use to study partial coverages
4. are we efficiently sampling both moves (i think so) and exchanges
5. Add forces to qe_out_info()
6. Finish distance rule by applying max displacement
7. learn about umbrella sampling
8. preferential sampling - surface/subsurface parameter (we kind of already do this with distance rules)
9. check behavior at 0 K and mu = 0 eV
10. show that high temperatures "heal" structures (e.g. move O from STO surface elsewhere and watch it find its way home) or correctly predict molecular transformations (e.g. 2 H2O2 <-> 2 H2O + O2)
11. improve sampling efficiency
12. replace moving steps with relax steps
13. force-bias or smart MC

list
1. grow nanoparticles (ti and o)
2. verify code is doing correct thing - set mu's in region where the bulk is stable, start with SrO and TiO2 surfaces and see if they grow more layers
3. select points at the edges and corners on sto phase diagram and grow surface
4. will ferroelectric (pto) grow different surfaces on each side
5. clean up code
6. implement simulated annealing
7. add forces to .axsf file
8. NP with high chemical potentials
9. add a distance rule to atom adding for surface simulations (min_dist = ri(A) + ri(B) and max_dist = value)
10. add a coordination number rule for NPs
11. SrO and TiO2 surfaces of SrTiO3 with different mu conditions (high Sr/low rest, high Ti/low rest, high O/low rest, high Sr & O/low Ti, high Ti & O/low Sr, high Sr & Ti/low O, high all)
12. make mu charts for each element
13. single atom relaxation for each adding step
14. mopac for initial scf calculations

readme
1. all_run.sh : runscript
2. el_list.txt : list of unique element indices, symbols, and masses
3. structure.txt : .xsf formatted file containing initial structure
4. /plot : contains plotting codes
