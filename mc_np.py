from mc import mc

class mc_np (mc):
	def __init__(self, T = 0, pace = 0, xsf = None):
		mc.__init__(self, T, pace, xsf) 
		
	# uvt update structure for nano particles
	# overrides the parent method
    def uvt_new_structure(self,xsf,el,act_p,bvo) : # el is of el_info class # act_p defines probability for different actions, [0]: move, [1]: swap, [2]: jump, [3]: add, [4]: remove
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