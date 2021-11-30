"""this module defines monte carlo objects and operations """
import random
import copy
import numpy as np
#from io import xsf_info
from structure import Element_info
#from bv import bv

from structure import Structure

# global variable definitions
ry_ev = 13.605693009
kb = 8.6173303e-5 # ev / k

class MonteCarlo :
    """class for representing monte carlo operations"""
    def __init__(self, T = 0, pace = 0, xsf = None) :
        self.T              = T                            # System temperature 
        self.T_max          = T                            # Max system temperature
        self.pace           = pace                       # Max displacement 
        self.pace_max       = 0.3                        # Max displacement limit    
        self.nvt_run_cnt    = 1                          # Canonical ensemble, running count, originally 0
        self.uvt_run_cnt    = 1                          # Grand canonical ensemble, ensemble running count, originally 0
        self.check_acc      = 25                          # Max displacement update rate, originally 0
        self.nvt_acc        = 0                          # Canonical ensemble, number of acceptance
        self.at_mv_num      = 0                          # Number of atoms to move
        self.at_mv_ind      = np.array([]).astype('int') # Indices of atoms to move
        self.at_mv_vec      = np.zeros((0, 3))           # Displacement of atoms
        if xsf is None:
            self.current_xsf        = Structure()                 # xsf in current iteration
            self.proposed_xsf        = Structure()                 # xsf in proposed iteration
            self.opt_xsf        = Structure()                 # xsf associated with the lowest energy
        else:
            self.current_xsf     = xsf.copy()
            self.proposed_xsf     = xsf.copy()
            self.opt_xsf     = xsf.copy()
        self.curr_en        = 100.0                        # Energy in current iteration
        self.curr_g         = 100.0                        # Free energy in current iteration
        self.opt_en         = 100.0                        # Lowest energy
        self.opt_g          = 100.0                        # Lowest free energy
        self.uvt_at_rm      = 0                          # Grand canonical ensemble, index of the atom to be removed
        self.uvt_at_add      = 0                          # Grand canonical ensemble, index of the atom to be added
        self.uvt_act        = 0                          # Grand canonical ensemble, 0: move, 1: swap, 2: jump, 3: add, 4: remove
        self.uvt_at_exc_num = 0                          # Grand canonical ensemble, number of exchanged atoms, can only be 0, 1, -1
        self.uvt_el_exc     = 0                          # Grand canonical ensemble, index of the element to be exchanged
        
    def copy(self) :
        return copy.deepcopy(self)
    
    # propose new coords
    def propose_coords(self) :
        # number of atoms to move
        self.at_mv_num = np.random.randint(1, self.current_xsf.atom_num + 1)
        if self.at_mv_num == self.current_xsf.atom_num :
            self.at_mv_ind = np.array(range(self.at_mv_num))
        else :
             self.at_mv_ind = np.array(random.sample(range(self.current_xsf.atom_num), self.at_mv_num))
        self.at_mv_vec = (np.random.rand(self.at_mv_num, 3) * 2 - 1) * self.pace

        self.proposed_xsf = self.current_xsf.copy()
        for i, ind in enumerate(self.at_mv_ind) :
            self.proposed_xsf.atom_coords[ind, :] += self.at_mv_vec[i, :]
    
    # update pace
    def update_pace(self) :
        if float(self.nvt_acc) / self.check_acc < 0.2 :
            self.pace /= 2.0
        elif float(self.nvt_acc) / self.check_acc > 0.4 :
            self.pace *= 1.618

        if self.pace > self.pace_max :
            self.pace = self.pace_max
        self.nvt_acc = 0
    
    def uvt_propose_structure(self, el, act_p, bvo) : # el is of el_info class # act_p defines probability for different actions, [0]: move, [1]: swap, [2]: jump, [3]: add, [4]: remove      
        '''
            uvt proposed structure based on random moves, deletes, inserts
        '''
        act_pp = np.array(act_p)
        # shorthand for the current structure, and also to be compatible with how Rob wrote the code
        # before I functionalized it and made it more like OOP
        xsf = self.current_xsf
        #------------------------------------------adjust act_p-----------------------------------------------------------------
        # avoid swapping action if only one element is removable(swappable)
        if act_pp[1] > 0 :
            el_swap_num = 0
            for i in range(len(self.current_xsf.atoms_swap)) :
                if len(self.current_xsf.atoms_swap[i]) > 0 :
                    el_swap_num += 1
            if el_swap_num <= 1 :
                act_pp[1] = 0
        # avoid jumping action if no appropriate site or no appropriate atom
        if act_pp[2] > 0 :
            # choose atom to jump
            at_neighbor_list = bvo.at_all_nn(xsf)
            at_neighbor_pref = np.zeros(self.current_xsf.atom_num).astype('int')
            for i in range(self.current.xsf.atom_num) :
                el_ind = self.current_xsf.at_type[i]
                at_neighbor_pref[i] = el.pref_nn[el_ind]
            weight = np.power((at_neighbor_list - at_neighbor_pref),4).astype('float')
            if np.sum(weight) != 0 :
                weight /= np.sum(weight)
                jump_at_ind = np.random.choice(range(xsf.atom_num), 1, p=weight)[0]
                jump_el_ind = self.current.xsf.at_type[jump_at_ind]
                """choose site to jump to"""
                for i in np.arange(xsf.vol * 1000) : 
                    jump_vec = np.zeros(3)
                    jump_vec += np.random.rand() * xsf.lat_vecs[0]
                    jump_vec += np.random.rand() * xsf.lat_vecs[1]
                    jump_vec += (np.random.rand() * (xsf.c_max - xsf.c_min) + xsf.c_min) / np.linalg.norm(xsf.lat_vecs[2]) * xsf.lat_vecs[2]
                    jump_neighbor = bvo.at_single_nn(xsf, jump_at_ind, jump_vec)
                    if (jump_neighbor == el.pref_nn[jump_el_ind]) : 
                        break
                if i >= xsf.vol * 1000 - 1 :
                    act_pp[2] = 0
            else :
                act_pp[2] = 0
        # avoid adding if no elements are able to add to the system
        if act_pp[3] > 0 :
            if np.sum(el.p_add) == 0 :
                act_pp[3] = 0
        # also avoid adding if could not find a site based on coordination rule
        if act_pp[3] > 0 :
            at_add = xsf.atom_num # index of atom to be added
            #print("el.num:", np.arange(el.num))
            #print("el.p_add:", el.p_add)
            self.uvt_el_exc  = np.random.choice(el.num, 1, p=el.p_add)[0]   # find the element index
            dis = 0
            trial = 0
            while dis < el.r_min[self.uvt_el_exc] or dis > el.r_max[self.uvt_el_exc] : # control atom distance
                at_add_coord  = np.zeros(3)
                at_add_coord += np.random.rand() * xsf.lat_vecs[0] 
                at_add_coord += np.random.rand() * xsf.lat_vecs[1]
                at_add_coord += (np.random.rand() * (xsf.c_max - xsf.c_min) + xsf.c_min) / np.linalg.norm(xsf.lat_vecs[2]) * xsf.lat_vecs[2]
                dis = min([np.linalg.norm(at_add_coord - xsf.atom_coords[ind]) for ind in range(xsf.atom_num)])
                trial += 1
                if trial >= 100000 :
                    break
            # don't add an atom if the distance to the nearest atom is too small or too large
            # this prevents the expected atom from very adding into the headspace far above the surface
            if trial >= 100000 and (dis < el.r_min[self.uvt_el_exc] or dis > el.r_max[self.uvt_el_exc]) :
                act_pp[3] = 0
        # avoid removing action if no removable atoms
        if act_pp[4] > 0 :
            if len(xsf.atoms_removable) == 0 : 
                act_pp[4] = 0
        # normalize act_p, and make it accumalate probability
        act_pp = act_pp / float(np.sum(act_pp))
        for i in range(len(act_pp) - 1) :
            act_pp[i + 1] += act_pp[i]
        #----------------------------------------end adjust act_p--------------------------------------------------------------
        # generate a random number between 0 and 1
        cndt = np.random.rand()

        #-------------move atoms--------------
        if cndt < act_pp[0] :    
            self.uvt_act        = 0
            self.uvt_at_exc_num = 0
            self.uvt_el_exc     = 0
            self.propose_coords()
            
        #-------------swap atoms-------------
        elif cndt < act_pp[1] : 
            self.uvt_act        = 1
            self.uvt_at_exc_num = 0
            # element to be swapped
            swap_el_1 = np.random.randint(el.num)
            while len(xsf.atoms_swap[swap_el_1]) == 0 :
                swap_el_1 = np.random.randint(el.num)
            swap_el_2 = np.random.randint(el.num)
            while swap_el_2 == swap_el_1 or len(xsf.atoms_swap[swap_el_2]) == 0 :
                swap_el_2 = np.random.randint(el.num)
            # atom to be swapped
            swap_at_1 = random.choice(current_xsf.atoms_swap[swap_el_1])
            swap_at_2 = random.choice(current_xsf.atoms_swap[swap_el_2])
            self.proposed_xsf.atom_coords[swap_at_1, :] = np.array(self.current_xsf.atom_coords[swap_at_2, :])
            self.proposed_xsf.atom_coords[swap_at_2, :] = np.array(self.current_xsf.atom_coords[swap_at_1, :])
            
        #-------------make atom jump, apply coord rule------------
        elif cndt < act_pp[2]:
            self.uvt_act        = 2
            self.uvt_at_exc_num = 0
            self.proposed_xsf.atom_coords[jump_at_ind, :] = np.array(jump_vec)
            
        #--------------------add one atom--------------------------------
        ### Need to move site searching part into adjusting p part, why??
        elif cndt < act_pp[3] : 
            self.uvt_act        = 3
            self.uvt_at_exc_num = 1
            self.proposed_xsf.atom_coords = np.vstack((self.current_xsf.atom_coords, at_add_coord))   # add coordinates to xsf
            self.proposed_xsf.atom_type  = np.append(self.current_xsf.atom_type, self.uvt_el_exc)  # add atom to the atom list
            self.proposed_xsf.atoms_rm   = np.append(self.current_xsf.atoms_removable, at_add)             # add atom to removable atom array
            self.proposed_xsf.atoms_swap[self.uvt_el_exc].append(at_add)              # add atom to swappable atom list
            self.proposed_xsf.num_each_element[self.uvt_el_exc] += 1                   # increase the number of that element
            self.proposed_xsf.atom_num += 1                                         # increase the total number of atoms
            
        #-----------------------remove one atom----------------------------
        else :             
            self.uvt_act        = 4
            self.uvt_at_exc_num = -1
            at_rm = random.choice(self.current_xsf.atoms_removable)                           # index of atom to be removed
            self.uvt_el_exc = self.current_xsf.atom_type[at_rm]                        # element index
            self.proposed_xsf.atom_coords = np.delete(xsf.atom_coords, at_rm, 0)   # remove the coordinates
            self.proposed_xsf.atom_type = np.delete(xsf.atom_type, at_rm, 0)     # remove the atom from atoms list
            self.proposed_xsf.num_each_element[self.uvt_el_exc] -= 1              # decrease the number of that element
            self.proposed_xsf.atom_num -= 1                                    # decrease the total number of atoms
            self.proposed_xsf.atoms_removable = np.append(xsf.atoms_removable[xsf.atoms_removable < at_rm], xsf.atoms_removable[xsf.atoms_removable > at_rm] - 1) # remove the atom from removable atoms
            # I don't think we always require that the removed atom
            # be in the swappable list, so this throws a ValueError if x is not in list
            if at_rm in self.proposed_xsf.atoms_swap[self.uvt_el_exc]:
                self.proposed_xsf.atoms_swap[self.uvt_el_exc].remove(at_rm)         # remove the atom from swappable atoms
            for i in range(len(self.proposed_xsf.atoms_swap)) :                        # remove the atom from swappable atoms
                for j in range(len(self.proposed_xsf.atoms_swap[i])) :
                    if self.proposed_xsf.atoms_swap[i][j] > at_rm :
                        self.proposed_xsf.atoms_swap[i][j] -= 1
        return self.proposed_xsf
  
    # Calculate free energy and probability threshhold of accepting steps
    # en in unit of eV
    def get_free_g_p(self, en, el, mu_list) : 
        free_g_new = en
        # calculate free energy
        for i in range(el.num) :
            free_g_new -= mu_list[i] * self.proposed_xsf.num_each_element[i]
        # update thermal de broglie wavelengths
        el.update_therm_db(self.T)
        # get thermal de broglie wavelengths of exchange element
        exc_therm_db = el.therm_db[self.uvt_el_exc]
        exc_el_num = self.proposed_xsf.num_each_element[self.uvt_el_exc] - (self.uvt_at_exc_num - 1) / 2
        if self.uvt_act <= 4 :
            exp_coef = np.exp(-(free_g_new - self.curr_g) / (self.T * kb))
        else : 
            print('Wrong action number for get_free_g_p! Undefined action!')
            exit()
        prob_acc = np.minimum(1, exp_coef * (self.proposed_xsf.vol / exc_therm_db**3 / exc_el_num )**self.uvt_at_exc_num)
        return free_g_new, prob_acc

    def update_T_const(self, iter, period) :
        self.T = self.T_max

    def update_T_linear(self, iter, period) : 
        self.T = self.T_max + float(1 - self.T_max) / period * (iter % period)

    def update_T_exp(self, iter, period) : 
        self.T = self.T-max * (1 / float(self.T_max))**( float(iter % period) / (period-1))

    def update_T_quadratic(self, iter, period) : 
        self.T = float(self.T_max - 1) / (period-1)**2 * (period - iter%period - 1)**2 + 1

    # canonical acceptance condition, return 1 if accepted, 0 otherwise;  energy in units eV
    def nvt_mc(self, en) :
        self.nvt_run_cnt += 1
        if self.nvt_run_cnt % self.check_acc == 0 :
                self.update_pace()
        if en < self.curr_en or np.random.uniform() < np.exp(-(en - self.curr_en) / (self.T * kb)) :
            self.curr_en = en
            self.nvt_acc += 1
            if en < self.opt_en :
                self.opt_en = en
                self.opt_xsf = self.proposed_xsf.copy()
            return 1
        else :       
            return 0
            
    # Grand canonical acceptance condition, return 1 if accepted, 0 otherwise
    def uvt_mc(self, en, el, mu_list) :
        self.uvt_run_cnt += 1
        rand = np.random.rand()
        free_g, prob_acc = self.get_free_g_p(en, el, mu_list)
        if self.uvt_act == 0 : # move
            self.nvt_run_cnt += 1
        elif self.uvt_act > 4 : 
            print('Wrong action number! Undefined action!')
            exit()
            if self.nvt_run_cnt % self.check_acc == 0:
                self.update_pace()
            if rand <= prob_acc: # accept
                if self.uvt_act == 0: # move
                    self.nvt_acc += 1
                self.curr_g = free_g
                self.current_xsf = self.proposed_xsf.copy()
                if free_g < self.opt_g : 
                    self.opt_g = free_g
                    self.opt_xsf = self.proposed_xsf.copy()
                return 1
            else:
                return 0
    # Some thing might need to do: instead of fixing the action probability, separate acc into acc_a, acc_r, acc_m, and adjust action probability based on these numbers, (basicly acc_a and acc_r), i.e. if acc_a is large, meaning that adding atom is favored, we should increase the action probability of add atoms
    
#  class that is a child of an MC process that adds the functionality to load proposed structures
#  from a lammps dump file
class MonteCarloLammps( MonteCarlo ):
    def __init__(self, T = 0, pace = 0, xsf = None):
        MonteCarlo.__init__(self, T, pace, xsf)
        
    def uvt_propose_structure_lammps(self, el, lmp, dump_file, step_max):
        """
            Use lammps object to load a structure from a set of dump files
            This structure will have been accepted as part of a markov
            chain with a classical hamiltonian
            It should be a close approximation of the QM-accurate states
        """
        #step_interval = 10.0
        step = np.random.randint(step_max)
        #step = step*step_interval
        print("step:", step)
        command = "read_dump " + dump_file + " " + str(step) + " x y z purge yes replace no add yes box no scaled no format xyz"
        print("LAMMPS COMMAND:", command)
        lmp.command(command)
        natoms = lmp.get_natoms()
        print("natoms", natoms)
        self.proposed_xsf = Structure(filename = None, el = el, natoms = natoms)
        # this assumes that the lattice vectors to do change from one snapshot to the next
        # i.e. no variable cells
        self.proposed_xsf.lat_vecs = self.current_xsf.lat_vecs.copy()
        self.proposed_xsf.atom_coords = lmp.numpy.extract_atom("x")
        self.proposed_xsf.atom_type = lmp.numpy.extract_atom("type").flatten()
        print("Atom coords:", self.proposed_xsf.atom_coords)
        print("Atom types:", self.proposed_xsf.atom_type)
        # print("One atom's coord:", 
        #    self.proposed_xsf.atom_coords[1][0], self.proposed_xsf.atom_coords[1][1], self.proposed_xsf.atom_coords[1][2])
        # print("Atom types:", self.proposed_xsf.atom_type)
        
        lmp.command("group groupSi type 1")
        lmp.command("group groupC type 2")
        lmp.command("variable nSi equal count(groupSi)")
        lmp.command("variable nC equal count(groupC)")
        
        np.append( self.proposed_xsf.num_each_element, lmp.numpy.extract_variable("nSi") )
        np.append( self.proposed_xsf.num_each_element, lmp.numpy.extract_variable("nC") )
        return self.proposed_xsf
        
def main_test():
    '''
    main function for debugging lammps functionality
    '''
    from lammps import lammps
    el_filename = "el_list.txt"
    lmp_init   = "gcmc.init"
    dump_file      = "gcmc.dump"
    step_max       = 25
    T = 500
    niter    = 100
    mc_run = MonteCarloLammps(T, 0.5)
    el = Element_info(el_filename, T) 
    energy_list = []
    lmp = lammps()
    lmp.file(lmp_init)
    for i in range(niter) :
        mc_run.uvt_propose_structure_lammps(el, lmp, dump_file, step_max)
        lmp.command("run 0")
        new_en = lmp.extract_compute("thermo_pe",0,0)
        print("Energy:", new_en)
        energy_list.append(new_en)
    print("energies::", energy_list)
        
# if running this as main, create a instance of MC for testing
if __name__ == "__main__":
    main_test()
