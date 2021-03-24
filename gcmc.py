#!/root/miniconda3/bin/python -i
# preceding line should have path for Python on your machine

# gcmc.py
# Purpose: mimic operation GCMC via Python
# Syntax:  gcmc.py in.gcmc
#          in.gcmc = LAMMPS input script

from __future__ import print_function
import sys,random,math

# set these parameters
# make sure neigh skin (in in.gcmc) > 2*deltamove

kb = 8.6173303e-5			# ev / k
random.seed(27848)

class FixGCMC : 
    def __init__(self, infile, lmp, fixId, type = 1, nmoves = 1, nexchanges = 2, gcmc_region = "all", group = "all"):
        self.deltamove = 0.1
        self.T  = 1000
        self.kT = kb*self.T
        
        self.ntotal = nmoves + nexchanges
        act_p = [nmoves/ntotal, 0.5*nexchanges/ntotal, 0.5*nexchanges/ntotal]
        act_cdf = []
        sum = 0
        for p in act_p:
            sum += p
            act_cdf.append( sum )
        print("cumulative action: ", act_cdf)
        self.gcmc_region = gcmc_region
        self.group = group
        self.seed = seed
               
        self.lmp = lmp        
        self.tmp = lammps()
        self.tmp.file(infile)
        # the temp is going to store the coordinates which are minimized
        # if the minimization isn't worth it do nothing
        # else replace the lmp.x array with the temp one

        self.tmp.command("variable etol equal 1.0e-8")
        self.tmp.command("variable ftol equal 1.0e-8")
        self.tmp.command("variable maxiter equal 500")
        self.tmp.command("variable maxeval equal 500")
        
        self.ngas_var = "ngas_" + group
        # e.g. 'variable nSi equal n(groupSi)'
        self.tmp.command("variable " + self.ngas_var + " equal count(" + group + ")")
        self.tmp.command("variable e equal pe")
        # we want to overwrite the dump and log to have seperate ones for temp and lmp
        tmp.command("dump tmpdump all xyz 10 " + fixId + "_tmp.xyz")
        tmp.command("log " + fixId + "_tmp.lammps")
        
        self.lmp.command("variable e equal pe")

        # run 0 to get energy of perfect lattice
        # emin = minimum energy

        self.lmp.command("run 0")

        natoms = lmp.extract_global("natoms",0)
        self.emin = lmp.extract_compute("thermo_pe",0,0) / natoms
        self.lmp.command("variable emin equal $e")
        
        lmp.command("variable elast equal $e")
        lmp.command("thermo_style custom step v_emin v_elast pe")
        lmp.command("run 0")
        x = lmp.extract_atom("x",3)
        lmp.command("variable elast equal $e")
        
        tmp.command("variable elast equal $e")
                  
        estart = lmp.extract_compute("thermo_pe",0,0) / natoms
     
    def run(self):
        # loop over Monte Carlo moves
        # extract x after every run, in case reneighboring changed ptr in LAMMPS
        self.elast = self.estart
        self.naccept = [0,0,0]
        was_accepted = False
        for i in range(self.ntotal):
          self.set_positions(self.lmp, self.tmp)
          self.update_ngas()
          cndt = random.random()
          if cndt <= self.act_cdf[0]:
            k = 0
            was_accepted =  self.move_attempt(self.tmp):
          elif cndt <= self.act_cdf[1]:
            k = 1
            was_accepted = self.insert_attempt(self.tmp)
          else:
            k = 2
            was_accepted = self.delete_attempt(self.tmp)
            
          if was_accepted:
            self.elast = self.e
            self.set_positions(self.tmp, self.lmp)
            self.lmp.command("variable elast equal $e")
            self.naccept[k] += 1         
               
    def set_positions(lmpFrom, lmpTo):
        '''
            Get the lmpTo system udpated from the lmpFrom system
            by updating position and velocity
            Normally this is set_positions(self.lmp, self.tmp)
            to set the temporary simulation box
        '''
        x = lmpFrom.gather_atoms("x", 1, 3)
        lmpTo.scatter_atoms("x", 1, 3, x)
        
        v = lmpFrom.gathejr_atoms("v", 1, 3)
        lmpTo.scatter_atoms("v", 1, 3, v)       

    def move_attempt(self, tmp):
        natoms = tmp.extract_global("natoms",0)
        iatom = random.randrange(0,natoms)
        x0 = x[iatom][0]
        y0 = x[iatom][1]
        z0 = x[iatom][2]

        x[iatom][0] += self.deltamove * (2*random.random()-1)
        x[iatom][1] += self.deltamove * (2*random.random()-1)
        x[iatom][2] += self.deltamove * (2*random.random()-1)
        
        tmp.command("run 1 pre no post no")
        x = tmp.extract_atom("x",3)
        self.e = tmp.extract_compute("thermo_pe",0,0)

        if self.e <= self.elast:
          return True
        elif random.random() <= math.exp((elast-e)/kT):
          return True
        else:
          x[iatom][0] = x0
          x[iatom][1] = y0
          x[iatom][2] = z0
          return False
        
    def insert_attempt(self, tmp):
        cndt = random.random()
        tmp.command("create_atoms " + str(self.type) + " random 1 " + str(self.lmp_seed) + " " + self.gcmc_region)
        
        tmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")        
        tmp.command("run 1 pre no post no")
        
        x = tmp.gather_atoms("x",1,3)
        self.e = tmp.extract_compute("thermo_pe",0,0)
        
        if self.e <= self.elast or random.random() <= math.exp((self.elast-self.e)/self.kT)*(self.zz*self.volume)/(self.ngas+1):
          return True
        else:
          return False
          
    def delete_attempt(self, tmp):
        natoms = self.lmp.extract_global("natoms",0)
        iatom = random.randrange(0, natoms)
        tmp.command("group exclude delete")
        tmp.command("group exclude id " + str(iatom))
        tmp.command("delete_atoms group exclude")
        
        tmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")        
        tmp.command("run 1 pre no post no")
        
        x = tmp.gather_atoms("x",1,3)
        self.e = tmp.extract_compute("thermo_pe",0,0)
        
        if self.e <= self.elast or random.random() <= self.ngas*math.exp((self.elast-self.e)/self.kT)/(self.zz*self.volume):
          return True
        else:
          return False
          
    def update_ngas(self, tmp):
        '''
            ngas is the number of atoms in the group
            for GCMC
        '''
        self.ngas = tmp.extract_variable(self.ngas_var)        
        
def final_stat(self, gcmc):
    # final energy and stats
    gcmc.lmp.command("variable nbuild equal nbuild")
    nbuild = gcmc.lmp.extract_variable("nbuild",None,0)
    natoms = gcmc.lmp.extract_global("natoms",0)

    gcmc.lmp.command("run 0")
    estop = gcmc.lmp.extract_compute("thermo_pe",0,0) / natoms

    print("MC stats:")
    print("  starting energy =", gcmc.estart)
    print("  final energy =", estop)
    #print("  minimum energy of perfect lattice =", gcmc.emin)
    print("  accepted MC moves =", gcmc.naccept)
    print("  neighbor list rebuilds =", nbuild)
        
def disorder_sys(self, lmp):
    # disorder the system
    # estart = initial energy
    self.deltaperturb = 0.5
    x = lmp.extract_atom("x",3)
    for i in range(natoms):
        x[i][0] += deltaperturb * (2*random.random()-1)
        x[i][1] += deltaperturb * (2*random.random()-1)
        x[i][2] += deltaperturb * (2*random.random()-1)
        
def main():
    # parse command line
    argv = sys.argv
    if len(argv) != 2:
      print("Syntax: gcmc.py in.gcmc")
      sys.exit()
    infile = sys.argv[1]
    from lammps import lammps
    lmp = lammps()
    # just sets up MC problem
    # run infile all at once
    lmp.file(infile)
    
    gcmcSi = FixGCMC(infile, lmp, "gcmcSi", type = 1)
    gcmcC = FixGCMC(infile, lmp, "gcmcC", type = 2)
    for n in range(1000):
        gcmc1.run()
        gcmc2.run()
        lmp.command("run 1")
    final_stats(gcmc)
    
if __name__ == "__main__":
    main()