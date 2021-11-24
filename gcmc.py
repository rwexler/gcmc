#!/root/miniconda3/bin/python -i
# preceding line should have path for Python on your machine

# gcmc.py
# Purpose: mimic operation GCMC via Python
# Syntax:  gcmc.py in.gcmc
#          in.gcmc = LAMMPS input script

from __future__ import print_function
import sys,random,math

from ctypes import c_double
from lammps import lammps
import math
# set these parameters
# make sure neigh skin (in in.gcmc) > 2*deltamove

kb = 8.6173303e-5			# ev / k
PI = 3.1415
HPLANCK = 4.135667e-3
#MVV2E = 
# conversion factor between mvv and eV

random.seed(27848)

class FixGCMC : 
    def __init__(self, infile, lmp, fixId, type = 1, nmoves = 1, nexchanges = 2, gcmc_region = "NULL", group = "all"):
        self.deltamove = 0.5
        self.T  = 1000
        self.kT = kb*self.T
        
        self.ntotal = nmoves + nexchanges
        act_p = [nmoves/self.ntotal, 0.5*nexchanges/self.ntotal, 0.5*nexchanges/self.ntotal]
        self.act_cdf = []
        sum = 0
        for p in act_p:
            sum += p
            self.act_cdf.append( sum )
        print("cumulative action: ", self.act_cdf)
        self.type = type
        self.gcmc_region = gcmc_region
        self.group = group
        self.seed = 12345

        self.beta = 1.0/self.kT
        #self.lambda = sqrt(HPLANCK*HPLANCK/(2.0*PI*self.mass*MVV2E*kb*self.T));
        if type == 1:
            chemName = "Si"
            chem_name = "si"
        else:
            chemName = "C"
            chem_name = "c"
        self.chemical_potential = lmp.extract_variable(chem_name + "_mu")
        self.mass = lmp.extract_variable("mass"+chemName)
        self.wavelength = lmp.extract_variable("broglie"+chemName)
        self.zz = math.exp(self.beta*self.chemical_potential)/(math.pow(self.wavelength,3.0))
           

        self.naccept = [0,0,0]
        self.nattempted = [0.1,0.1,0.1]
               
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
        
        # we want to overwrite the dump and log to have seperate ones for temp and lmp
        self.tmp.command("dump tmpdump all xyz 10 " + fixId + "_tmp.xyz")
        self.tmp.command("log " + fixId + "_tmp.lammps")
        
        self.tmp.command("variable e equal pe")
        self.lmp.command("variable e equal pe")

        # run 0 to get energy of perfect lattice
        # emin = minimum energy
        self.tmp.command("run 0")
        self.lmp.command("run 0")

        self.volume = 11931.703 #lmp.extract_compute("vol", 0, 0) #self.get_volume()

        natoms = self.tmp.extract_global("natoms",0)
        self.e = self.tmp.extract_compute("thermo_pe",0,0)
        #self.lmp.command("variable emin equal $e")
        
        self.tmp.command("variable elast equal $e")
        self.tmp.command("variable estart equal $e")
        self.tmp.command("thermo_style custom step v_estart v_elast pe")
        self.tmp.command("thermo_modify lost warn")
        self.tmp.command("run 0")

        self.lmp.command("variable elast equal $e")
        self.lmp.command("variable estart equal $e")
        self.lmp.command("thermo_style custom step v_estart v_elast pe")
        self.lmp.command("run 0")
        #x = lmp.extract_atom("x",3)
        #self.lmp.command("variable elast equal $e")
        
        self.tmp.command("variable elast equal $e")
        self.lmp.command("variable elast equal $e")

        self.estart = self.tmp.extract_compute("thermo_pe",0,0)
        self.elast = self.estart

        self.update_ngas(self.tmp)
     
    def run(self):
        # loop over Monte Carlo moves
        # extract x after every run, in case reneighboring changed ptr in LAMMPS
        was_accepted = False
        for i in range(self.ntotal):
          self.set_positions(self.lmp, self.tmp)
          self.update_ngas(self.tmp)
          cndt = random.random()
          if cndt <= self.act_cdf[0]:
            k = 0
            was_accepted =  self.move_attempt(self.tmp)
          elif cndt <= self.act_cdf[1]:
            k = 1
            was_accepted = self.insert_attempt(self.tmp)
          else:
            k = 2
            was_accepted = self.delete_attempt(self.tmp)
            
          if was_accepted:
            print("ACCEPTED")
            self.elast = self.e
            self.set_positions(self.tmp, self.lmp)
            #self.lmp.command("variable elast equal $e")
            self.naccept[k] += 1         
          self.nattempted[k] += 1
               
    def set_positions(self, lmpFrom, lmpTo):
        '''
            Get the lmpTo system udpated from the lmpFrom system
            by updating position and velocity
            Normally this is set_positions(self.lmp, self.tmp)
            to set the temporary simulation box
        '''
        x = lmpFrom.gather_atoms("x", 1, 3)
        lmpTo.scatter_atoms("x", 1, 3, x)
        
        v = lmpFrom.gather_atoms("v", 1, 3)
        lmpTo.scatter_atoms("v", 1, 3, v)
        # make sure energies are updated
        lmpTo.command("run 0 post no")

    def move_attempt(self, tmp):
        print("ATTEMPT MC MOVE")
        natoms = tmp.extract_global("natoms",0)
        #x = (natoms*3*c_double)( tmp.gather_atoms("x", 1, 3) )
        x = tmp.extract_atom("x")
        iatom = random.randrange(0,natoms)
        #print("x[iatom][0]: ", x[iatom][0])
        x0 = x[iatom][0]
        y0 = x[iatom][1]
        z0 = x[iatom][2]

        x[iatom][0] += self.deltamove * (2*random.random()-1)
        x[iatom][1] += self.deltamove * (2*random.random()-1)
        x[iatom][2] += self.deltamove * (2*random.random()-1)
        
        tmp.command("run 0 pre no post no")
        #x = tmp.extract_atom("x",3)
        self.e = tmp.extract_compute("thermo_pe",0,0)
        print("self.e: ", self.e)
        print("self.elast: ", self.elast)
        if self.e <= self.elast:
          return True
        elif random.random() <= math.exp((self.elast-self.e)/self.kT):
          return True
        else:
          x[iatom][0] = x0
          x[iatom][1] = y0
          x[iatom][2] = z0
          return False
        
    def insert_attempt(self, tmp):
        print("ATTEMPT INSERT")
        cndt = random.random()
        current_natoms = tmp.get_natoms()

        tmp.command("create_atoms " + str(self.type) + " random 1 " + str(self.seed) + " " + self.gcmc_region)
        try:
            tmp.command("group ADD delete")
        except:
            print("no group ADD")
        tmp.command("group ADD id " + str(current_natoms))
        tmp.command("delete_atoms overlap 1.0 ADD all")        
        # the extra atom was deleted
        if tmp.get_natoms() == current_natoms:
            return False
        tmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")        
        tmp.command("run 0 pre no post no")
        
        #x = tmp.gather_atoms("x",1,3)
        self.e = tmp.extract_compute("thermo_pe",0,0)
        
        if self.e <= self.elast or random.random() <= math.exp((self.elast-self.e)/self.kT)*(self.zz*self.volume)/(self.ngas+1):
          return True
        else:
          return False
          
    def delete_attempt(self, tmp):
        print("ATTEMPT DELETE")
        natoms = self.lmp.extract_global("natoms",0)
        iatom = random.randrange(0, natoms)
        # the group exclude may not exist yet
        try:    
            tmp.command("group exclude delete")
        except:
            print("no group exclude")
        #tmp.command("group clear exclude")
        tmp.command("group exclude id " + str(iatom))
        tmp.command("delete_atoms group exclude")
        
        tmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")        
        tmp.command("run 0 pre no post no")
        
        #x = tmp.gather_atoms("x",1,3)
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
        
def final_stats(gcmc):
    # final energy and stats
    gcmc.lmp.command("variable nbuild equal nbuild")
    nbuild = gcmc.lmp.extract_variable("nbuild",None,0)
    #natoms = gcmc.lmp.extract_global("natoms",0)

    gcmc.lmp.command("run 0")
    estop = gcmc.lmp.extract_compute("thermo_pe",0,0)
    paccept = [gcmc.naccept[k]/gcmc.nattempted[k] for k in range(3)]

    print("MC stats:")
    print("  starting energy =", gcmc.estart)
    print("  final energy =", estop)
    #print("  minimum energy of perfect lattice =", gcmc.emin)
    print("  accepted MC moves =", gcmc.naccept)
    print("  percent of MC moves accepted = ", paccept)
    print("  neighbor list rebuilds =", nbuild)
        
def disorder_sys(lmp):
    # disorder the system
    # estart = initial energy
    deltaperturb = 0.5
    x = lmp.extract_atom("x",3)
    for i in range(lmp.extract_global("natoms", 0)):
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
    lmp = lammps()
    # just sets up MC problem
    # run infile all at once
    lmp.file(infile)

    #disorder_sys(lmp)
    
    gcmcSi = FixGCMC(infile, lmp, "gcmcSi", type = 1, nexchanges = 2)
    gcmcC = FixGCMC(infile, lmp, "gcmcC", type = 2, nexchanges = 2)
    for n in range(1000):
        gcmcSi.run()
        gcmcC.run()
        lmp.command("run 1 pre no post no")
    final_stats(gcmcSi)
    
if __name__ == "__main__":
    main()
