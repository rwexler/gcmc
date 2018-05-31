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
import subprocess

####################
# READ INPUT FILES #
####################
xsf_filename   = sys.argv[1] # read xsf filename
el_filename    = sys.argv[2] # read element list filename

#############################
# SET SIMULATION PARAMETERS #
#############################
niter    = 1000
max_disp = 0.05                               # angstroms
T_move   = 500                                # kelvin
ry_ev    = 13.605693009
bohr_ang = 0.52917721067
buf_len  = 3.5                                # length above surface within which atoms can be added
mu_list  = [-989.926, -428.350]               # ag, o - ag/ag2o
#mu_list  = [-989.866, -428.468]               # ag, o - ag-rich
#mu_list  = [-989.985, -428.231]               # ag, o - o-rich
act_p    = np.array([1e-5, 1e-5, 1e-5, 1, 1]) # probablity of taking different actions
					      # [0]: move, [1]: swap, [2]: jump, [3]: add, [4]: remove
fail_en  = 999.
nproc    = 72
nkdiv    = 2

###########################
# GET ELEMENT INFORMATION #
###########################
el = el_info() 
el.pop_attr(el_filename,T_move)

##############################
# GET STRUCTURAL INFORMATION #
##############################
xsf = xsf_info()
xsf.pop_attr(xsf_filename, el, buf_len)

####################################
# GET NEAREST NEIGHBOR INFORMATION #
####################################
bvo = bv()
bvo.init(xsf, el)

###################################
# INSTANTIATE MONTE CARLO ROUTINE #
###################################
mc_run = mc()
mc_run.init(T_move, max_disp, xsf)

###################################
# RUN GRAND CANONICAL MONTE CARLO #
###################################
os.system('mkdir -p temp')                                              # make temp directory for qe calculations
os.chdir('temp')                                                        # enter temp
log_file              = init_log('log.dat')                             # initialize log file
axsf_opt_file         = init_axsf('coord_opt.axsf', niter, xsf)         # "        " axsf file recording optimized structure
axsf_new_file         = init_axsf('coord_new.axsf', niter, xsf)         # "                            " structure created in current iteration
axsf_accept_file      = init_axsf('coord_accept.axsf', niter, xsf)      # "                                      " accepted in current iteration
axsf_failed_file      = init_axsf('coord_failed.axsf', niter, xsf)      # initialize axsf file recording structure failed in qe
axsf_failed_iter_file = init_axsf('coord_failed_iter.axsf', niter, xsf) # initialize axsf file recording structure failed in qe
failed_cnt = 0
for i in range(niter) :
	# attempt uvt action and store xsf attributes in xsf_new
	if i == 0 : 
		# alway start with move
		xsf = mc_run.uvt_new_structure(xsf, el, np.array([1,0,0,0,0]), bvo) 
	else :
		xsf = mc_run.uvt_new_structure(xsf, el, act_p, bvo) 

	# make input file
	make_qe_in('qe.in', xsf, el)

	# if not move step, then replace scf with relax
	if mc_run.uvt_act != 0 :
		os.system('sed -i "s/scf/relax/g" qe.in')
	
	# if number of atoms is smaller than or equal to 6, make sure it has 36 bands
	# note that in the future it should be evaluated based on number of electrons
	# and the threshold for the number of bands should be -ndiag in qe
	if xsf.at_num <= 6 :
		os.system('sed -i "/&SYSTEM/a nbnd = 36" qe.in')
	
	# calculate and get total energy
	call_qe = 'mpiexec.hydra -np ' + str(nproc) + ' ../bin/pw.x -nk ' + str(nkdiv) + ' -i qe.in > qe.out'
	subprocess.call(call_qe, shell = True)
	qe_out = qe_out_info('qe.out')

	# get energy and forces from qe
	if os.popen('grep ! qe.out').read() == '' :
		# qe failed at first scf step
		new_en = fail_en
		mc_run.new_xsf.at_force = np.zeros((mc_run.new_xsf.at_num, 3))
	else :
		# get energy from qe
		new_en = qe_out.get_final_en()
		# get forces from qe
		mc_run.new_xsf.at_force = qe_out.get_forces(mc_run.new_xsf.at_num) * ry_ev / bohr_ang # convert forces from ry/bohr to ev/ang

	# if not move step and qe does not fail at first step, update atomic coordinates
	if ( os.popen('grep ! qe.out').read() != '' and mc_run.uvt_act != 0 ) :
		mc_run.new_xsf.at_coord = qe_out.get_coord(mc_run.new_xsf.at_num)

	# update T
	mc_run.update_T_const(i - failed_cnt, 3000)

	# decide whether or not to accept uvt action 
	# note that old_xsf is changed to new_xsf if accepted
	accept = mc_run.uvt_mc(new_en, el, mu_list)

	# calculate free energy 
	free_en, _ = mc_run.get_free_g_p(new_en, el, mu_list)

	# if step not accepted, copy attributes from old (previous) xsf to xsf
	if accept == 0 :
		xsf = mc_run.old_xsf.copy()
	# othewise copy new xsf to xsf
	else :
		xsf = mc_run.new_xsf.copy()

	# update logs if no qe error
	if new_en != fail_en :
		# write energies, number of accepted steps, and acceptance rate to log file
		upd_log(log_file, i - failed_cnt, free_en, mc_run)
		# write atomic coordinates to axsf file
		upd_axsf(axsf_opt_file, i - failed_cnt, mc_run.opt_xsf, el)
		upd_axsf(axsf_new_file, i - failed_cnt, mc_run.new_xsf, el)
		upd_axsf(axsf_accept_file, i - failed_cnt, mc_run.old_xsf, el)
	else :
		failed_cnt += 1
		upd_axsf(axsf_failed_file, failed_cnt, mc_run.new_xsf, el)
		upd_axsf(axsf_failed_iter_file, i, mc_run.new_xsf, el)

log_file.close()
axsf_opt_file.close()
axsf_new_file.close()
axsf_accept_file.close()
axsf_failed_file.close()
axsf_failed_iter_file.close()
os.chdir('../')
