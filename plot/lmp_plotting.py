import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lmp_utils import *

def plot_comp(df, filename, makePDF = True):
    '''
    Plot the composition of the Markov Chain against timestep
    '''
    # what is the difference between df['nSi'] and df.loc['nSi']
    nSi = df['nSi'].values 
    nC = df['nC'].values 
    indices = df['indx'].values 
    
    fig, ax = plt.subplots(figsize=(8, 4))
    #x = np.arange(0, num_structs, 1)
    ax.plot(indices, nSi, 'k--', linewidth=1.5, label='Si')
    ax.plot(indices, nC, 'r--', linewidth=1.5, label='C')
    ax.legend(loc='right')
    ax.set_xlabel('Step Number')
    ax.set_ylabel('Composition')
    ax.set_title('SiC GCMC')
    #plt.show()
    plt.savefig( generate_plot_filename(filename, "_comp") )

def plot_stats(df, filename, makePDF = True):
    '''
    Plot the other stats from lammps log file against timestep
    Potential Energy, Reduced Potential, Acceptance Pecentage or Absolute
    '''
    #new_en = np.asarray(new_en)
    #en_curr = np.asarray(en_curr)
    #en_low = np.asarray(en_low)
    #acc = np.asarray(acc)
    #pace = np.asarray(pace)
    
    # define a cutoff for which the system is still equilibrating
    # cutoff = 50
    
    poteng = df['PE'].values #[cutoff:] 
    redPot = df['redPot'].values
    nacc_insSi = df['nacc_insSi'].values
    nacc_delSi = df['nacc_delSi'].values
    nacc_insC = df['nacc_insC'].values
    nacc_delC = df['nacc_delC'].values
    indices = df['indx'].values
    
    nSi = df['nSi'].values 
    nC = df['nC'].values 
    comp = nSi/nC
    
    # plot contents of log file
    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8, 4))

    # plot energies
    ax0.plot(indices, poteng, 'b--', label = 'Potential')
    #ax0.legend(loc="upper right")
    ax0.set_xlabel('Step Number')
    ax0.set_ylabel('Potential Energy (eV)', color = 'b')
    ax0.tick_params('y', labelcolor='b')
    #en_min = np.min(poteng)
    #en_max = np.max(poteng)
    #ax0.set_ylim(en_min, en_max)
    ax3 = ax0.twinx()
    ax3.plot(indices, redPot, 'r--', label = 'Reduced Potential')
    ax3.set_ylabel('Grand Potential Energy (eV)', color = 'r')
    ax3.tick_params('y', labelcolor='r')
    
    ax1.plot(indices, comp, 'k--', linewidth=1.5, label='Si/C')
    ax1.legend(loc='right')
    ax1.set_xlabel('Step Number')
    ax1.set_ylabel('% Composition')
    ax1.set_title('Ratio of Composition nSi/nC')

    fig.tight_layout()
    plt.savefig( generate_plot_filename(filename, "_stats") )

def plot_accs(df, filename, makePDF = True):
    nacc_insSi = df['nacc_insSi'].values 
    nacc_delSi = df['nacc_delSi'].values
    nacc_insC = df['nacc_insC'].values
    nacc_delC = df['nacc_delC'].values
    indices = df['indx'].values
    
    prcnt_insSi = nacc_insSi/indices
    prcnt_delSi = nacc_delSi/indices
    prcnt_insC = nacc_insC/indices
    prcnt_delC = nacc_delC/indices
    
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 4))
    # plot number of accepted steps and the percentage of accepted steps
    ax1.plot(indices, nacc_insSi, 'k-', label = "Inserted Si")
    ax1.plot(indices, nacc_delSi, 'b-', label = "Deleted Si")
    ax1.plot(indices, nacc_insC, 'r-', label = "Inserted C")
    ax1.plot(indices, nacc_delC, 'm-', label = "Deleted C")   
    ax1.set_xlabel('Step Number')
    ax1.set_ylabel('Number of Accepted Steps', color = 'k')
    ax1.legend(loc="lower right")

    ax2.plot(indices, prcnt_insSi, 'k--', label = "Inserted Si")
    ax2.plot(indices, prcnt_delSi, 'b--', label = "Deleted Si")
    ax2.plot(indices, prcnt_insC, 'r--', label = "Inserted C")
    ax2.plot(indices, prcnt_delC, 'm--', label = "Deleted C")
    ax2.set_xlabel('Step Number')
    ax2.set_ylabel('Percent of Steps Accepted', color = 'k')
    ax2.legend(loc="upper right")

    fig.tight_layout()
    plt.savefig( generate_plot_filename(filename, "_accs") )
    
def plot_phase_diagram(df, filename, makePDF = True):    
    dmu_line = np.linspace(-2, 0, 101)
    pref_phase = np.zeros((dmuOs.shape[0], ))

    for index, row in df.iterrows() :
        # calculate gammaO
        gammaO = row.loc['GammaO']
        # calculate phi
        phi = row.loc['Phi']        
        # get surface area
        A = row.loc['A']
        omega = np.zeros((dmu_line.shape[0], ))
        for i, dmu in enumerate(dmu_line) :            
            # calculate omega in J/m^2
            omega[i] = (1. / (2. * A)) * (phi + gammaO * dmu) * 16.0218
            # 1 eV/A^2 = 16.0218 J/m^2
        # check if preferred phase
        if index == 0 :
            omega_min = omega
        else :
            for i in range(dmu_line.shape[0]) :
                if omega[i] < omega_min[i] :
                    pref_phase[i] = index
                    omega_min[i] = omega[i]

        # plot surface energy line
        plt.plot(dmu_line, omega, alpha = 0.6, zorder = 1)
        # plot first surface energy line
        if index == 0 :
            plt.plot(dmu_line, omega, 'b', lw = 2, zorder = 2)

    # plot minimum surface energy and print preferred phases
    plt.plot(dmu_line, omega_min, 'k', lw = 2, zorder = 2)
    for i in range(dmu_line.shape[0]) :
        print(dmu_line[i], pref_phase[i], omega_min[i])

    # plot bulk stability region
    min = -3
    max = 1
    xmin = np.repeat(min, 10)
    xmax = np.repeat(max, 10)
    y = np.linspace(0, 5, 10) 
    plt.plot(xmin, y, 'k:', zorder = 3)
    plt.plot(xmax, y, 'k:', zorder = 3)
    # plot format
    plt.xlim(-2, 0)
    plt.ylim(0.8, 1.1)
    plt.xlabel(r'$\Delta \mu_{\rm O}$ (eV)')
    plt.ylabel(r'Surface energy (J/m$^2$)')
    #plt.show()
    plt.savefig( generate_plot_filename(filename, "_phase_diagram") )

def plot_spd_stats(df, filename, makePDF = True):
    f, (ax1, ax2) = plt.subplots(1, 2)

    ratio = df['nSi'].values / df['nC'].values
    init_ratio = ratio[0]
    ax1.hist(ratio, zorder = 1)
    ax1.axvline(init_ratio, color = 'k', ls = '--', zorder = 2)
    ax1.set_xlabel(r'$n_{\rm Ag} / n_{\rm O}$')
    ax1.set_yticks([])

    df['Phi'].plot(kind = 'density')
    ax2.set_xlabel(r'Surface energy (J/m$^2$) at $\Delta \mu_{\rm O} = 0$ eV')
    ax2.set_yticks([])
    ax2.set_ylabel('')
    #plt.show()
    plt.savefig( generate_plot_filename(filename, "_spd_stats") )