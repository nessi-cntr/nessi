import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from ReadCNTR import write_input_file
#----------------------------------------------------------------------
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z
#----------------------------------------------------------------------
def GenSysParams(Nsites,HoppingT,HubbardU,MuChem,Beta,FermBos=-1):
    sysparams = {'Nsites': Nsites,
                 'HoppingT': HoppingT,
                 'HubbardU': HubbardU,
                 'MuChem': MuChem,
                 'beta': Beta,
                 'FermBos': FermBos
        }
    return sysparams
#----------------------------------------------------------------------
def GenRampParams(RampSite,RampW0):
    rampparams = {'RampSite': RampSite,
                  'RampW0': RampW0
        }
    return rampparams
#----------------------------------------------------------------------
def GenSolverParams(Nt,Ntau,dt,MatsMaxIter=100,MatsMaxErr=1.0e-8,BootstrapMaxIter=20,
                        BootstrapMaxErr=1.0e-8,CorrectorSteps=3,SolverOrder=5):

    solverparams = {'Nt': Nt,
                    'Ntau': Ntau,
                    'dt': dt,
                    'MatsMaxIter': MatsMaxIter,
                    'MatsMaxErr': MatsMaxErr,
                    'BootstrapMaxIter': BootstrapMaxIter,
                    'BootstrapMaxErr': BootstrapMaxErr,
                    'CorrectorSteps': CorrectorSteps,
                    'SolverOrder': SolverOrder
                    }
    return solverparams
#----------------------------------------------------------------------
def GenInputFile(input_file,sysparams,rampparams,solverparams):
    #inp = {**sysparams, **rampparams, **solverparams}
    inp = merge_two_dicts(sysparams, rampparams)
    inp = merge_two_dicts(inp, solverparams)
    write_input_file(inp, input_file)
#----------------------------------------------------------------------
def RunHubbardChain(sysparams,rampparams,solverparams,approx_sigma,output_file,
                        input_file='',log_file='',runpath='./'):
    ns = sysparams['Nsites']
    U = sysparams['HubbardU']
    Nt = solverparams['Nt']
    Ntau = solverparams['Ntau']
    dt = solverparams['dt']
    if len(input_file) == 0:
        flin = './inp/hubbardchain_N{}_U{}_Nt{}_Ntau{}_dt{}.inp'.format(ns,U,Nt,Ntau,dt)
    else:
        flin = input_file
    GenInputFile(flin,sysparams,rampparams,solverparams)
    if approx_sigma.lower() == '2b':
        prog = runpath + '/exe/hubbard_chain_2b.x'
    elif approx_sigma.lower() == 'gw':
        prog = runpath + '/exe/hubbard_chain_gw.x'
    elif approx_sigma.lower() == 'tpp':
        prog = runpath + '/exe/hubbard_chain_tpp.x'
    else:
        print('[ERROR] Self-energy approximation not recognized.')
        exit()

    log_flag = ''
    if len(log_file) > 0:
        log_flag = ' > ' + log_file

    os.system(prog + ' ' + flin + ' ' + output_file + log_flag)
#----------------------------------------------------------------------

if __name__ == '__main__':

    showref = False

    if len(sys.argv) > 1:
        flref = sys.argv[1]
        print('[Info] Using reference data from:' + flref)
        showref = True

    Nsites=2  # number of sites
    thop=1.0  # hopping amplitude
    U=2.0     # Hubbard repulsion
    mu=0.0    # chemical potential
    beta=20.0 # inverse temperature
    sysparams = GenSysParams(Nsites,thop,U,mu,beta)

    RampSite=1    # on which site the ramp is applied
    RampW0 = 1.0  # magnitude of ramp (w0)
    rampparams = GenRampParams(RampSite,RampW0)

    Nt = 200    # number of time steps
    Ntau = 400  # number of points on Matsubara branch
    dt = 0.025  # time step
    solverparams = GenSolverParams(Nt,Ntau,dt,BootstrapMaxErr=1.0e-10)

    fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)

    sigma_approx = ['2B','GW','TPP']

    color_palette = ['blue', 'orange', 'green']

    ax2.set_xlabel(r'$t$')
    ax2.set_ylabel('energy')
    ax1.set_ylabel(r'$n_1(t)$')

    ax1.set_xlim(0.0,Nt*dt)
    ax2.set_xlim(0.0,Nt*dt)

    for i,approx in enumerate(sigma_approx):
        output_file = 'out/hubbardchain_' + approx
        RunHubbardChain(sysparams,rampparams,solverparams,approx,output_file)

        data = np.loadtxt(output_file + '_occupation.dat')
        t = data[:,0]
        ax1.plot(t,data[:,1],label=approx,c=color_palette[i])

        t,Ekin,Epot,Etot = np.loadtxt(output_file + '_energy.dat',unpack=True)
        ax2.plot(t,Ekin,c=color_palette[i])
        ax2.plot(t,Etot,c=color_palette[i],linestyle='--')

    if showref:
        tr,occr,Ekinr = np.loadtxt(flref,unpack=True)
        ax1.plot(tr,occr,c='k',label='exact',linestyle=':')
        ax2.plot(tr,Ekinr,c='k',linestyle=':')

    ax1.legend(loc='upper right',frameon=False)
    plt.show()
