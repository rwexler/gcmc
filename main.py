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
from lammps import lammps

####################
# READ INPUT FILES #
####################
xsf_filename   = sys.argv[1] # read xsf filename
el_filename    = sys.argv[2] # read element list filename

lmp_init       = ""
dump_file      = ""
step_max       = 5000
biasProposals = False       # don't use lammps, just the basic random structure proposal

#############################
# SET SIMULATION PARAMETERS #
#############################
niter    = 1000
max_disp = 0.0                                # angstroms
T_move   = 500                                # kelvin
ry_ev    = 13.605693009
bohr_ang = 0.52917721067
buf_len  = 2.5                                # length above surface within which atoms can be added
#mu_list  = [-989.926, -428.156]               # ag, o - p/p0 = 1.e+0
#mu_list  = [-989.926, -428.451]               # ag, o - p/p0 = 1.e-10
#mu_list  = [-989.926, -428.333]               # ag, o - p/p0 = 1.e-6
#mu_list  = [-989.926, -428.215]               # ag, o - p/p0 = 1.e-2
mu_list  = [-989.926, -428.096]               # ag, o - p/p0 = 1.e+2
act_p    = np.array([0.3, 0, 0, 0.35, 0.35]) # probablity of taking different actions
					      # [0]: move, [1]: swap, [2]: jump, [3]: add, [4]: remove
fail_en  = 999
nproc    = 144
nkdiv    = 1
ndiag    = 144

useQE = False						# parameter to check if Quantum Espresso i.e. DFT calculations should be used

###########################
# GET ELEMENT INFORMATION #
###########################
el = Element_info(el_filename, T_move) 

##############################
# GET STRUCTURAL INFORMATION #
##############################
if useQE:
	xsf = Structure_xsf(xsf_filename, el, buf_len)
else
	xsf = None
####################################
# GET NEAREST NEIGHBOR INFORMATION #
####################################
bvo = BondValence()

###################################
# INSTANTIATE MONTE CARLO ROUTINE #
###################################
mc_run = MonteCarlo(T_move, max_disp, xsf)

lmp = lammps()
lmp.file(lmp_init)
		
###################################
# RUN GRAND CANONICAL MONTE CARLO #
###################################
os.system('mkdir -p tmp')                                              # make temporary directory for qe calculations
os.chdir('tmp')                                                        # enter temp
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
		mc_run.uvt_propose_structure(el, np.array([1,0,0,0,0]), bvo) 
	else :
		if biasProposals:
            # load biased proposals from LAMMPS dump files, proposals that are part of a markov chain with the classical hamiltonian
            mc_run.uvt_propose_structure_lammps(el, lmp, dump_file, step_max)
        else:
            mc_run.uvt_propose_structure(el, act_p, bvo) 

    # get energy and structure from QE and XSF
    # make input file
    make_qe_in('qe.in', mc_run.proposed_xsf, el)

    # if not move step, then replace scf with relax
    if mc_run.uvt_act != 0 :
        os.system('sed -i "s/scf/relax/g" qe.in')
		
    # if number of atoms is smaller than or equal to 6, make sure it has 36 bands
    # note that in the future it should be evaluated based on number of electrons
    # and the threshold for the number of bands should be -ndiag in qe
    if mc_run.proposed_xsf.atom_num <= 6 :
        os.system('sed -i "/&SYSTEM/a nbnd = 36" qe.in')
		
    # calculate and get total energy by calling QE
    call_qe = 'mpiexec_mpt -np ' + str(nproc) + ' ../bin/pw.x -nk ' + str(nkdiv) + ' -ndiag ' + str(ndiag) + ' -i qe.in > qe.out'
    subprocess.call(call_qe, shell = True)
    qe_out = qe_out_info('qe.out')

    # get energy and forces from qe
    if os.popen('grep ! qe.out').read() == '' :
        # qe failed at first scf step
        new_en = fail_en
        mc_run.proposed_xsf.atom_forces = np.zeros((mc_run.proposed_xsf.atom_num, 3))
    else :
        # get energy from qe
        new_en = qe_out.get_final_en()
        # get forces from qe
        mc_run.proposed_xsf.atom_forces = qe_out.get_forces(mc_run.proposed_xsf.atom_num) * ry_ev / bohr_ang # convert forces from ry/bohr to ev/ang

    # if not move step and qe does not fail at first step, update atomic coordinates
    if ( os.popen('grep ! qe.out').read() != '' and mc_run.uvt_act != 0 ) :
        mc_run.proposed_xsf.atom_coords = qe_out.get_coord(mc_run.proposed_xsf.atom_num)
		      
	# decide whether or not to accept uvt action 
	# note that current_xsf is changed to proposed_xsf if accepted if it was a move
	accept = mc_run.uvt_mc(new_en, el, mu_list)

	# calculate free energy 
	free_en, _ = mc_run.get_free_g_p(new_en, el, mu_list)

	# if step is accepted set current xsf to the proposed (and now accepted) structure
    if accept == 1 :
		mc_run.current_xsf = mc_run.proposed_xsf.copy()

	# update logs if no qe error
	if new_en != fail_en :
		# write energies, number of accepted steps, and acceptance rate to log file
		upd_log(log_file, i - failed_cnt, free_en, mc_run)
		# write atomic coordinates to axsf file
		upd_axsf(axsf_opt_file, i - failed_cnt, mc_run.opt_xsf, el)
		upd_axsf(axsf_new_file, i - failed_cnt, mc_run.proposed_xsf, el)
		upd_axsf(axsf_accept_file, i - failed_cnt, mc_run.current_xsf, el)
	else :
		failed_cnt += 1
		upd_axsf(axsf_failed_file, failed_cnt, mc_run.proposed_xsf, el)
		upd_axsf(axsf_failed_iter_file, i, mc_run.proposed_xsf, el)

log_file.close()
axsf_opt_file.close()
axsf_new_file.close()
axsf_accept_file.close()
axsf_failed_file.close()
axsf_failed_iter_file.close()
os.chdir('../')
