import sys
import copy
from io import xsf_info, el_info
from io import qe_out_info, make_qe_in
from io import init_log, upd_log
from io import init_axsf, upd_axsf
from mc import mc
from bv import bv
import os
import numpy as np

xsf_filename = sys.argv[1] # read xsf filename from command line
el_filename = sys.argv[2] # read element list filename

# set simulation parameters
niter = 1000
max_disp = 0.05 # angstroms
T_move = 1 # kelvin
ry_ev = 13.605693009
bohr_ang = 0.52917721067
buf_len = 2.0 # length above surface within which atoms can be added
mu_list = [0, 0, 0] # sr, ti, o
act_p = np.array([1,1,1,0,0]) # probablity of taking different actions, [0]: move, [1]: swap, [2]: jump, [3]: add, [4]: remove

# get element info
el = el_info() # instantiates el_info object
el.pop_attr(el_filename,T_move)

# get info from xsf file
xsf = xsf_info() # instantiate xsf_info objects
xsf.pop_attr(xsf_filename, el, buf_len) # populate attributes in xsf_info object

# create bv object
bvo = bv()
bvo.init(xsf)

# instantiate mc object
mc_run = mc()
mc_run.init(T_move, max_disp, xsf)

# run mc simulation
os.system('mkdir -p temp') # make temp directory for qe calculations
os.chdir('temp') # enter temp
log_file = init_log('log.dat') # initialize log file
axsf_opt_file    = init_axsf('coord_opt.axsf', niter, xsf)    # initialize axsf file recording optimized structure
axsf_curr_file   = init_axsf('coord_curr.axsf', niter, xsf)   # initialize axsf file recording structure created in current iteration
axsf_accept_file = init_axsf('coord_accept.axsf', niter, xsf) # initialize axsf file recording structure accepted in current iteration
i = 0
while i < niter :
	# attempt uvt action and store xsf attributes in xsf_new
	if i == 0 : 
		# alway start with moving
		xsf = mc_run.uvt_new_structure(xsf, el, np.array([1,0,0,0,0]), bvo) 
	else :
		xsf = mc_run.uvt_new_structure(xsf, el, act_p, bvo) 

	# make input file
	make_qe_in('qe.in', xsf, el)
	
	# calculate and get total energy
	if xsf.at_num <= 2 :
		os.system('mpiexec.hydra -np 4 ../bin/pw.x -i qe.in > qe.out') # execute qe
	elif xsf.at_num <= 6 :
		os.system('mpiexec.hydra -np 9 ../bin/pw.x -i qe.in > qe.out')
	else :
		os.system('mpiexec.hydra -np 36 ../bin/pw.x -i qe.in > qe.out')
	qe_out = qe_out_info('qe.out')
	# get energy from qe
	new_en = qe_out.get_final_en() * ry_ev # convert final energy from ry to ev
	# get forces from qe
	mc_run.new_xsf.at_force = qe_out.get_forces(mc_run.new_xsf.at_num) * ry_ev / bohr_ang # convert forces from ry/bohr to ev/ang
	xsf.at_force = np.array(mc_run.new_xsf.at_force)

	# update T
	mc_run.update_T_const(i, 3000)

	# decide whether or not to accept uvt action, 
	accept = mc_run.uvt_mc(new_en, el, mu_list)

	# calculate free energy 
	free_en, _ = mc_run.get_free_g_p(new_en, el, mu_list)

	# if step not accepted, copy attributes from old (previous) xsf to xsf
	if accept == 0 :
		xsf = mc_run.old_xsf.copy()

	# update logs if no errors in qe running
	if new_en != 0.0 :
		# write energies, number of accepted steps, and acceptance rate to log file
		upd_log(log_file, i, free_en, mc_run)

		# write atomic coordinates to axsf file
		upd_axsf(axsf_opt_file, i, mc_run.opt_xsf, el)
		upd_axsf(axsf_curr_file, i, mc_run.new_xsf, el)
		upd_axsf(axsf_accept_file, i, mc_run.old_xsf, el)
	else :
		i -= 1
	
	i += 1

log_file.close()
axsf_file.close()
os.chdir('../')
