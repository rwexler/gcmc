"""this module defines monte carlo objects and operations """
##### : Need to implement/correct immediately
####  : Need to implement/correct when switching to another function
###   : A very rare case may cause error
##    : Better implementation
import random
import copy
import numpy as np
from io import xsf_info
from io import el_info
from bv import bv

# global variable definitions
ry_ev = 13.605693009 # ev
kb    = 8.6173303e-5 # ev / k

class mc :
    """class for representing monte carlo operations"""

    def __init__(self) :
        self.T              = 0.0                        # system temperature 
        self.T_max          = 0.0                        # max system temperature
        self.pace           = 0.0                        # max displacement 
        self.pace_max       = 0.0                        # max displacement limit 
        self.nvt_run_cnt    = 0                          # canonical ensemble, running count, how many steps are executed
        self.uvt_run_cnt    = 0                          # grand canonical ensemble, running count
        self.check_acc      = 0                          # max displacement update rate
        self.nvt_acc        = 0                          # canonical ensemble, how many steps are accepted
        self.at_mv_num      = 0                          # number of atoms to move
        self.at_mv_ind      = np.array([]).astype('int') # indices of atoms to move
        self.at_mv_vec      = np.zeros((0, 3))           # displacement of atoms
        self.old_xsf        = xsf_info()                 # xsf in the previous iteration
        self.new_xsf        = xsf_info()                 # xsf in the current iteration
        self.opt_xsf        = xsf_info()                 # xsf with the lowest energy
        self.curr_en        = 0.0                        # energy in the current iteration
        self.curr_g         = 0.0                        # free energy in the current iteration
        self.opt_en         = 0.0                        # lowest energy
        self.opt_g          = 0.0                        # lowest free energy
        self.uvt_at_rm      = 0                          # grand canonical ensemble, index of the atom to be removed
        self.uvt_at_ad      = 0                          # grand canonical ensemble, index of the atom to be added
        self.uvt_act        = 0                          # grand canonical ensemble, 0: move, 1: swap, 2: jump, 3: add, 4: rem
        self.uvt_at_exc_num = 0                          # grand canonical ensemble, number of exch atoms, b/w -1 and 1 # RBW
        self.uvt_el_exc     = 0                          # grand canonical ensemble, index of the element to be exchanged

    def copy(self) :
        """perform explicit copy of mc"""
        cp_self = mc()
        cp_self.T              = self.T            
        cp_self.T_max          = self.T_max
        cp_self.pace           = self.pace          
        cp_self.pace_max       = self.pace_max      
        cp_self.nvt_run_cnt    = self.nvt_run_cnt   
        cp_self.uvt_run_cnt    = self.uvt_run_cnt   
        cp_self.check_acc      = self.check_acc     
        cp_self.nvt_acc        = self.nvt_acc       
        cp_self.at_mv_num      = self.at_mv_num
        cp_self.at_mv_ind      = np.array(self.at_mv_ind)     
        cp_self.at_mv_vec      = np.array(self.at_mv_vec)     
        cp_self.old_xsf        = self.old_xsf.copy()       
        cp_self.new_xsf        = self.new_xsf.copy()       
        cp_self.opt_xsf        = self.opt_xsf.copy()       
        cp_self.curr_en        = self.curr_en       
        cp_self.curr_g         = self.curr_g        
        cp_self.opt_en         = self.opt_en        
        cp_self.opt_g          = self.opt_g         
        cp_self.uvt_at_rm      = self.uvt_at_rm 
        cp_self.uvt_at_ad      = self.uvt_at_ad 
        cp_self.uvt_act        = self.uvt_act       
        cp_self.uvt_at_exc_num = self.uvt_at_exc_num   
        cp_self.uvt_el_exc     = self.uvt_el_exc
        return cp_self
    
    def init(self, T, pace, xsf) :
        """initialize mc"""
        self.T           = T
        self.T_max       = T
        self.pace        = pace
        self.pace_max    = 0.3
        self.nvt_run_cnt = 1
        self.uvt_run_cnt = 1
        self.check_acc   = 25
        self.old_xsf     = xsf.copy()
        self.new_xsf     = xsf.copy()
        self.opt_xsf     = xsf.copy()
        self.curr_en     = 99.
        self.curr_g      = 99.
        self.opt_en      = 99.
        self.opt_g       = 99.

    def rand_mv(self, xsf) :
        """move atoms"""

        self.at_mv_num = np.random.randint(1, xsf.at_num + 1)                           # number of atoms to move
        if self.at_mv_num == xsf.at_num :                                               # if all atoms
            self.at_mv_ind = np.array(range(self.at_mv_num))                            # indices are 1 to number of atoms
        else :                                                                          # if not all atoms
            self.at_mv_ind = np.array(random.sample(range(xsf.at_num), self.at_mv_num)) # rand select indices
        self.at_mv_vec = (np.random.rand(self.at_mv_num, 3) * 2 - 1) * self.pace        # gen disps

    def new_coords(self, xsf) :
        """update atomic coordinates"""
        self.rand_mv(xsf)                                         # gen disps
        self.old_xsf = xsf.copy()                                 # copy xsf
        self.new_xsf = xsf.copy()
        for i, ind in enumerate(self.at_mv_ind) :                 # loop over atoms to move
            self.new_xsf.at_coord[ind, :] += self.at_mv_vec[i, :] # apply disps
    
    # update pace
    def update_pace(self) :
        if float(self.nvt_acc) / self.check_acc < 0.2 :
            self.pace /= 2.0
        elif float(self.nvt_acc) / self.check_acc > 0.4 :
            self.pace *= 1.618

        if self.pace > self.pace_max :
            self.pace = self.pace_max
        self.nvt_acc = 0
    
    def uvt_new_structure(self,xsf,el,act_p,bvo) :
        """update structure
                el contains element info
                act_p defines prob for diff act
                        [0] : move
                        [1] : swap
                        [2] : jump
                        [3] : add
                        [4] : remove
        """

        self.old_xsf = xsf.copy()
        self.new_xsf = xsf.copy()
        act_pp = np.array(act_p)  # editable copy of act_p

        ###############################
        # ADJUST ACTION PROBABILITIES #
        ###############################

        # avoid swapping action if only one element is swappable
        if act_pp[1] > 0 :                     # if possible to swap
            el_swap_num = 0                    # number of swappable elem
            for i in range(len(xsf.at_swap)) : # loop over swappable elem
                if len(xsf.at_swap[i]) > 0 :   # if, for that elem, there are swappable atoms
                    el_swap_num += 1           # incr number of swappable elem
            if el_swap_num <= 1 :              # if there are no swappable elem
                act_pp[1] = 0                  # avoid swapping

        # avoid jumping action if there are no suitable sites or atoms
        if act_pp[2] > 0 :                                                                 # if possible to jump

            # choose atom to jump
            at_neighbor_list = bvo.at_all_nn(xsf)                                          # nn for each atom
            at_neighbor_pref = np.zeros(xsf.at_num).astype('int')                          # pref nn
            for i in range(xsf.at_num) :                                                   # loop over atoms
                el_ind = xsf.at_type[i]                                                    # elem index
                at_neighbor_pref[i] = el.pref_nn[el_ind]                                   # get pref nn
            weight = np.power((at_neighbor_list - at_neighbor_pref), 4).astype('float')    # penalty for under/over-coord
            if np.sum(weight) != 0 :                                                       # if any atoms are penalized
                weight /= np.sum(weight)                                                   # normalize penalty
                jump_at_ind = np.random.choice(range(xsf.at_num), 1, p = weight)[0]        # rand select atom w/ prob = penalty
                jump_el_ind = xsf.at_type[jump_at_ind]                                     # what elem is it

                # choose site to jump to
                for i in np.arange(xsf.vol * 1000) :                                       # samp vol w/ avg gr sp of vol/1000
                    jump_vec = np.zeros(3)                                                 # init pos of new site
                    jump_vec += np.random.rand() * xsf.lat_vec[0]                          # rand select x pos
                    jump_vec += np.random.rand() * xsf.lat_vec[1]                          # y pos
                    jump_vec += (np.random.rand() * (xsf.c_max - xsf.c_min) + xsf.c_min) \
                        / np.linalg.norm(xsf.lat_vec[2]) * xsf.lat_vec[2]                  # z pos
                    jump_neighbor = bvo.at_single_nn(xsf, jump_at_ind, jump_vec)           # nn @ new site
                    if (jump_neighbor == el.pref_nn[jump_el_ind]) :                        # if nn = pref found 
                        break                                                              # allow jumping

                if i >= xsf.vol * 1000 - 1 :                                               # else
                    act_pp[2] = 0                                                          # avoid jumping
            else :                                                                         # if no atoms are penalized
                act_pp[2] = 0                                                              # avoid jumping

        # avoid adding action if there are no elements that are able to added to the system
        if act_pp[3] > 0 :             # if possible to add
            if np.sum(el.p_add) == 0 : # if sum of add prob = 0
                act_pp[3] = 0          # avoid adding

        # avoid adding action if a site could not be found that satisfies the coordination rule
        if act_pp[3] > 0 :                                                                # if possible to add
            at_ad = xsf.at_num                                                            # index of new atom
            self.uvt_el_exc = np.random.choice(range(el.num), 1, p = el.p_add)[0]         # choose elem
            dis = 0                                                                       # init min dist b/w new atom and old
            trial = 0                                                                     # init trial count
            while dis < el.r_min[self.uvt_el_exc] or dis > el.r_max[self.uvt_el_exc] :    # if min dist not w/in rmin and rmax
                at_ad_coord  = np.zeros(3)                                                # init pos of new atom
                at_ad_coord += np.random.rand() * xsf.lat_vec[0]                          # rand select x pos # RBW
                at_ad_coord += np.random.rand() * xsf.lat_vec[1]                          # y pos # RBW
                at_ad_coord += (np.random.rand() * (xsf.c_max - xsf.c_min) + xsf.c_min) \
                    / np.linalg.norm(xsf.lat_vec[2]) * xsf.lat_vec[2]                     # z pos
                dis = min([np.linalg.norm(at_ad_coord - xsf.at_coord[ind]) \
                               for ind in range(xsf.at_num)])                             # calc min dist
                trial += 1                                                                # count trial
                if trial >= 100000 :                                                      # if while not closed in 100k trials
                    break                                                                 # exit
            if trial >= 100000 and (dis < el.r_min[self.uvt_el_exc] \
                                        or dis > el.r_max[self.uvt_el_exc]) :             # if there were 100k trials
                act_pp[3] = 0                                                             # avoid adding

        # avoid removing action if there are no removable atoms
        if act_pp[4] > 0 :            # if possible to remove
            if len(xsf.at_rmb) == 0 : # if there are no removable atoms
                act_pp[4] = 0         # avoid removing

        # norm act_pp and conv it to a cum prob
        act_pp = act_pp / float(np.sum(act_pp)) # norm to 1
        for i in range(len(act_pp) - 1) :       # loop through probs
            act_pp[i + 1] += act_pp[i]          # add the prev prob

        ###################
        # PERFORM ACTIONS #
        ###################

        cndt = np.random.rand() # gen rand number b/w 0 and 1

        # move ---------------------------------------------
        if cndt < act_pp[0] :          # if move
            self.uvt_act        = 0    # set index of action
            self.uvt_at_exc_num = 0    # set number of atoms to exch
            self.uvt_el_exc     = 0    # index of elem to be exch
            self.new_coords(xsf)       # move atoms
            return self.new_xsf.copy()

        # swap ---------------------------------------------
        elif cndt < act_pp[1] :                                                                 # if swap
            self.uvt_act        = 1                                                             # set index of action
            self.uvt_at_exc_num = 0                                                             # set number of atoms to exch

            # elements to be swapped
            swap_el_1 = np.random.randint(el.num)                                               # rand select elem index
            while len(xsf.at_swap[swap_el_1]) == 0 :                                            # if no swap atoms of that elem
                swap_el_1 = np.random.randint(el.num)                                           # rand select new elem index
            swap_el_2 = np.random.randint(el.num)                                               # rand select 2nd elem index
            while swap_el_2 == swap_el_1 or len(xsf.at_swap[swap_el_2]) == 0 :                  # if elems are same or no swap atoms of 2nd elem
                swap_el_2 = np.random.randint(el.num)                                           # rand select new 2nd elem index

            # atoms to be swapped
            swap_at_1 = random.choice(xsf.at_swap[swap_el_1])                                   # rand select atom
            swap_at_2 = random.choice(xsf.at_swap[swap_el_2])                                   # rand select 2nd atom
            self.new_xsf.at_coord[swap_at_1, :] = np.array(self.old_xsf.at_coord[swap_at_2, :]) # replace 2nd atom with 1st
            self.new_xsf.at_coord[swap_at_2, :] = np.array(self.old_xsf.at_coord[swap_at_1, :]) # replace 1st atom with 2nd
            return self.new_xsf.copy()

        # jump ---------------------------------------------
        elif cndt < act_pp[2]:                                         # if jump
            self.uvt_act        = 2                                    # set index of action
            self.uvt_at_exc_num = 0                                    # set number of atoms to exch
            self.new_xsf.at_coord[jump_at_ind, :] = np.array(jump_vec) # move atom to jump site
            return self.new_xsf.copy()

        # add ----------------------------------------------
        elif cndt < act_pp[3] :                                             # if add
            self.uvt_act        = 3                                         # set index of action
            self.uvt_at_exc_num = 1                                         # set number of atoms to exch
            self.new_xsf.at_coord = np.vstack((xsf.at_coord, at_ad_coord))  # add coords to xsf
            self.new_xsf.at_type  = np.append(xsf.at_type, self.uvt_el_exc) # add atom to atom list
            self.new_xsf.at_rmb   = np.append(xsf.at_rmb, at_ad)            # add atom to removable atom array
            self.new_xsf.at_swap[self.uvt_el_exc].append(at_ad)             # add atom to swappable atom list
            self.new_xsf.el_each_num[self.uvt_el_exc] += 1                  # incr number of that elem
            self.new_xsf.at_num += 1                                        # incr total number of atoms
            return self.new_xsf.copy()

        # remove -------------------------------------------
        else :                                                                  # if remove
            self.uvt_act        = 4                                             # set index of action
            self.uvt_at_exc_num = -1                                            # set number of atoms to exch
            at_rm = random.choice(xsf.at_rmb)                                   # index of atom to remove
            self.uvt_el_exc = xsf.at_type[at_rm]                                # elem index
            self.new_xsf.at_coord = np.delete(xsf.at_coord, at_rm, 0)           # remove coords
            self.new_xsf.at_type = np.delete(xsf.at_type, at_rm, 0)             # remove atom
            self.new_xsf.el_each_num[self.uvt_el_exc] -= 1                      # decr number of that elem
            self.new_xsf.at_num -= 1                                            # decr total number of atoms
            self.new_xsf.at_rmb = np.append(xsf.at_rmb[xsf.at_rmb < at_rm], 
                                            xsf.at_rmb[xsf.at_rmb > at_rm] - 1) # change indices for removable atoms
            self.new_xsf.at_swap[self.uvt_el_exc].remove(at_rm)                 # remove atom from swappable atoms
            for i in range(len(self.new_xsf.at_swap)) :                         # loop through elems
                for j in range(len(self.new_xsf.at_swap[i])) :                  # loop through atoms of that elem
                    if self.new_xsf.at_swap[i][j] > at_rm :                     # if atom index > index of atom to remove
                        self.new_xsf.at_swap[i][j] -= 1                         # reduce index of those atoms
            return self.new_xsf.copy()

    # uvt update structure for nano particles
    def uvt_new_structure_np(self,xsf,el,act_p,bvo) : # el is of el_info class # act_p defines probability for different actions, [0]: move, [1]: swap, [2]: jump, [3]: add, [4]: remove
        self.old_xsf = xsf.copy()
        self.new_xsf = xsf.copy()
        act_pp = np.array(act_p)
        #------------------------------------------adjust act_p-----------------------------------------------------------------
        # avoid swapping action if only one element is removable(swappable)
        if act_pp[1] > 0 :
            el_swap_num = 0
            for i in range(len(xsf.at_swap)) :
                if len(xsf.at_swap[i]) > 0 :
                    el_swap_num += 1
            if el_swap_num <= 1 :
                act_pp[1] = 0
        # avoid jumpping action if no appropriate site or no appropriate atom
        if act_pp[2] > 0 :
            # choose atom to jump
            at_neighbor_list = bvo.at_all_nn(xsf)
            at_neighbor_pref = np.zeros(xsf.at_num).astype('int')
            for i in range(xsf.at_num) :
                el_ind = xsf.at_type[i]
                at_neighbor_pref[i] = el.pref_nn[el_ind]
            weight = np.power((at_neighbor_list - at_neighbor_pref),4).astype('float')
            if np.sum(weight) != 0 :
                weight /= np.sum(weight)
                jump_at_ind = np.random.choice(range(xsf.at_num), 1, p=weight)[0]
                jump_el_ind = xsf.at_type[jump_at_ind]
                """choose site to jump to"""
                for i in np.arange(xsf.vol * 1000) : 
                    jump_vec = np.zeros(3)
                    jump_vec += np.random.rand() * xsf.lat_vec[0]
                    jump_vec += np.random.rand() * xsf.lat_vec[1]
                    jump_vec += (np.random.rand() * (xsf.c_max - xsf.c_min) + xsf.c_min) / np.linalg.norm(xsf.lat_vec[2]) * xsf.lat_vec[2]
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
            at_ad = xsf.at_num # index of atom to be added
            self.uvt_el_exc  = np.random.choice(range(el.num), 1, p=el.p_add)[0]   # find the element index
            dis = 0
            trial = 0
            while dis < el.r_min[self.uvt_el_exc] or dis > el.r_max[self.uvt_el_exc] : # control atom distance
                at_ad_coord  = np.random.rand(3) * 2 - 1
                while np.linalg.norm(at_ad_coord) > 1 or np.linalg.norm(at_ad_coord) < xsf.r_min/xsf.r_max :
                    at_ad_coord  = np.random.rand(3) * 2 - 1
                at_ad_coord *= xsf.r_max
                dis = min([np.linalg.norm(at_ad_coord - xsf.at_coord[ind]) for ind in range(xsf.at_num)])
                trial += 1
                if trial >= 100000 :
                    break
            if trial >= 100000 and (dis < el.r_min[self.uvt_el_exc] or dis > el.r_max[self.uvt_el_exc]) :
                act_pp[3] = 0
        # avoid removing action if no removable atoms
        if act_pp[4] > 0 :
            if len(xsf.at_rmb) == 0 : 
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
            self.new_coords(xsf)
            self.rand_mv(xsf)
            return self.new_xsf.copy()

        #-------------swap atoms-------------
        elif cndt < act_pp[1] : 
            self.uvt_act        = 1
            self.uvt_at_exc_num = 0
            # element to be swapped
            swap_el_1 = np.random.randint(el.num)
            while len(xsf.at_swap[swap_el_1]) == 0 :
                swap_el_1 = np.random.randint(el.num)
            swap_el_2 = np.random.randint(el.num)
            while swap_el_2 == swap_el_1 or len(xsf.at_swap[swap_el_2]) == 0 :
                swap_el_2 = np.random.randint(el.num)
            # atom to be swapped
            swap_at_1 = random.choice(xsf.at_swap[swap_el_1])
            swap_at_2 = random.choice(xsf.at_swap[swap_el_2])
            self.new_xsf.at_coord[swap_at_1, :] = np.array(self.old_xsf.at_coord[swap_at_2, :])
            self.new_xsf.at_coord[swap_at_2, :] = np.array(self.old_xsf.at_coord[swap_at_1, :])
            return self.new_xsf.copy()

        #-------------make atom jump, apply coord rule------------
        elif cndt < act_pp[2]:
            self.uvt_act        = 2
            self.uvt_at_exc_num = 0
            self.new_xsf.at_coord[jump_at_ind, :] = np.array(jump_vec)
            return self.new_xsf.copy()

        #--------------------add one atom--------------------------------
        ### Need to move site searching part into adjusting p part
        elif cndt < act_pp[3] : 
            self.uvt_act        = 3
            self.uvt_at_exc_num = 1
            self.new_xsf.at_coord = np.vstack((xsf.at_coord, at_ad_coord))   # add coordinates to xsf
            self.new_xsf.at_type  = np.append(xsf.at_type, self.uvt_el_exc)  # add atom to the atom list
            self.new_xsf.at_rmb   = np.append(xsf.at_rmb, at_ad)             # add atom to removable atom array
            self.new_xsf.at_swap[self.uvt_el_exc].append(at_ad)              # add atom to swappable atom list
            self.new_xsf.el_each_num[self.uvt_el_exc] += 1                   # increase the number of that element
            self.new_xsf.at_num += 1                                         # increase the total number of atoms
            return self.new_xsf.copy()

        #-----------------------remove one atom----------------------------
        else :             
            self.uvt_act        = 4
            self.uvt_at_exc_num = -1
            at_rm = random.choice(xsf.at_rmb)                           # index of atom to be removed
            self.uvt_el_exc = xsf.at_type[at_rm]                        # element index
            self.new_xsf.at_coord = np.delete(xsf.at_coord, at_rm, 0)   # remove the coordinates
            self.new_xsf.at_type = np.delete(xsf.at_type, at_rm, 0)     # remove the atom from atoms list
            self.new_xsf.el_each_num[self.uvt_el_exc] -= 1              # decrease the number of that element
            self.new_xsf.at_num -= 1                                    # decrease the total number of atoms
            self.new_xsf.at_rmb = np.append(xsf.at_rmb[xsf.at_rmb < at_rm], xsf.at_rmb[xsf.at_rmb > at_rm] - 1) # remove the atom from removable atoms
            self.new_xsf.at_swap[self.uvt_el_exc].remove(at_rm)         # remove the atom from swappable atoms
            for i in range(len(self.new_xsf.at_swap)) :                        # remove the atom from swappable atoms
                for j in range(len(self.new_xsf.at_swap[i])) :
                    if self.new_xsf.at_swap[i][j] > at_rm :
                        self.new_xsf.at_swap[i][j] -= 1
            return self.new_xsf.copy()


    # canonical acceptance condition, return 1 if accepted, 0 otherwise
    # en in unit of eV
    def nvt_mc(self, en) :
        self.nvt_run_cnt += 1
        if en < self.curr_en :
            self.curr_en = en
            self.nvt_acc += 1
            if en < self.opt_en :
                self.opt_en = en
                self.opt_xsf = self.new_xsf.copy()
            if self.nvt_run_cnt % self.check_acc == 0 :
                self.update_pace()
            return 1
        elif np.random.uniform() < np.exp(-(en - self.curr_en) / (self.T * kb)) :
            self.curr_en = en
            self.nvt_acc += 1
            if self.nvt_run_cnt % self.check_acc == 0 :
                self.update_pace()
            return 1
        else :
            if self.nvt_run_cnt % self.check_acc == 0 :
                self.update_pace()
            return 0

    def get_free_g_p(self, en, el, mu_list) :
        """calc free energy (in ev) and acc prob"""
        free_g_new = en                                                    # init free energy as dft total energy
        for i in range(el.num) :                                           # loop through elem
            free_g_new -= mu_list[i] * self.new_xsf.el_each_num[i]         # subtract off chem pots for each elem
        el.update_therm_db(self.T)                                         # upd therm de broglie wavelengths
        exc_therm_db = el.therm_db[self.uvt_el_exc]                        # get therm de broglie wavelength of elem to be exch
        exc_el_num = self.new_xsf.el_each_num[self.uvt_el_exc] \
            - (self.uvt_at_exc_num - 1) / 2                                # equals n + 1 for add and n for remove
        if self.uvt_act <= 4 :                                             # if action exists
            exp_coef = np.exp(-(free_g_new - self.curr_g) / (self.T * kb)) # calc exp prefactor for acc prob
        else :                                                             # if action does not exist
            print 'error in get_free_g_p: action does not exist!'          # print error message
            exit()                                                         # exit
        prob_acc = np.minimum(1, exp_coef * (self.new_xsf.vol / exc_therm_db ** 3 / exc_el_num) ** self.uvt_at_exc_num)
        return free_g_new, prob_acc

    # adjust temp
    def update_T_const(self, iter, period) :
        """maintain const temp"""
        self.T = self.T_max

    def update_T_linear(self, iter, period) :
        """sawtooth temp prof"""
        self.T = self.T_max + float(1 - self.T_max) / period * (iter % period)

    def update_T_exp(self, iter, period) : 
        """exp sawtooth temp prof"""
        self.T = self.T_max * (1 / float(self.T_max)) ** ( float(iter % period) / (period - 1))

    def update_T_quadratic(self, iter, period) : 
        """quad sawtooth temp prof"""
        self.T = float(self.T_max - 1) / (period - 1) ** 2 * (period - iter % period - 1) ** 2 + 1

    def uvt_mc(self, en, el, mu_list) :
        """grand canonical acceptance condition
                accepted = 1
                rejected = 0
        """
        self.uvt_run_cnt += 1                                 # count run
        rand = np.random.rand()                               # select rand number b/w 0 and 1
        free_g, prob_acc = self.get_free_g_p(en, el, mu_list) # calculate g and acc prob

        if self.uvt_act == 0 :                                # if move
            self.nvt_run_cnt += 1                             # count run
            if rand <= prob_acc :                             # accept
                self.nvt_acc += 1                             # count acc step
                self.curr_g = free_g                          # update free energy
                self.old_xsf = self.new_xsf.copy()            # update structure
                if free_g < self.opt_g :                      # if lowest free energy
                    self.opt_g = free_g                       # update min free energy
                    self.opt_xsf = self.new_xsf.copy()        # update opt structure
                if self.nvt_run_cnt % self.check_acc == 0 :   # if appropriate
                    self.update_pace()                        # update max disp
                return 1
            else :                                            # reject
                if self.nvt_run_cnt % self.check_acc == 0 :   # if appropriate
                    self.update_pace()                        # update max disp
                return 0

        elif self.uvt_act == 1 or self.uvt_act == 2 :         # if swap or jump
            if rand <= prob_acc :                             # accept
                self.curr_g = free_g                          # update free energy
                self.old_xsf = self.new_xsf.copy()            # update structure
                if free_g < self.opt_g :                      # if lowest free energy
                    self.opt_g = free_g                       # update min free energy
                    self.opt_xsf = self.new_xsf.copy()        # update opt structure
                return 1
            else :                                            # reject
                return 0

        elif self.uvt_act == 3 or self.uvt_act == 4 :         # if exchange
            if rand <= prob_acc :                             # accept
                self.curr_g = free_g                          # update free energy
                self.old_xsf = self.new_xsf.copy()            # update structure
                if free_g < self.opt_g :                      # if lowest free energy
                    self.opt_g = free_g                       # update min free energy
                    self.opt_xsf = self.new_xsf.copy()        # update opt structure
                return 1
            else :                                            # reject
                return 0
        else :                                                # if action does not exist
            print 'error in uvt_mc: action does not exist!'   # print error message
            exit()                                            # exit
