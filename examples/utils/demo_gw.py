import sys
import os
import glob
import numpy as np
#import matplotlib.pyplot as plt
from ReadCNTR import write_input_file
#----------------------------------------------------------------------
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z
#----------------------------------------------------------------------
def GenSysParams(Nk,HoppingT,HubbardU,V,MuChem,Beta,FermBos=-1):
    sysparams = {'Nk': Nk,
                 'HoppingT': HoppingT,
                 'HubbardU': HubbardU,
                 'V': V,
                 'MuChem': MuChem,
                 'beta': Beta,
                 'FermBos': FermBos
        }
    return sysparams
#----------------------------------------------------------------------
def GenSolverParams(Nt,Ntau,dt,SaveGreen=False,MatsMaxIter=100,
                        MatsMaxErr= 1.0e-8,
                        BootstrapMaxIter= 20,BootstrapMaxErr=1.0e-8,
                        TimeMaxErr=1.0e-8,CorrectorSteps=10,SolverOrder=5):

    savegf = 0
    if SaveGreen:
        savegf = 1
    
    solverparams = {'Nt': Nt,
                    'Ntau': Ntau,
                    'dt': dt,
                    'MatsMaxIter': MatsMaxIter,
                    'MatsMaxErr': MatsMaxErr,
                    'BootstrapMaxIter': BootstrapMaxIter,
                    'BootstrapMaxErr': BootstrapMaxErr,
                    'TimeMaxErr': TimeMaxErr,
                    'CorrectorSteps': CorrectorSteps,
                    'SolverOrder': SolverOrder,
                    'SaveGreen': savegf
                        }
    return solverparams
#----------------------------------------------------------------------
def GenInputFile(input_file,sysparams,solverparams):
    inp = sysparams
    inp = merge_two_dicts(inp, solverparams)
    write_input_file(inp, input_file)
#----------------------------------------------------------------------
def Run(sysparams,solverparams,output_file,
                        input_file='',log_file='',runpath=''):
    ns = sysparams['Nk']
    U = sysparams['HubbardU']
    Nt = solverparams['Nt']
    Ntau = solverparams['Ntau']
    dt = solverparams['dt']
    if len(input_file) == 0:
        flin = './inp/hubbardchain_N{}_U{}_Nt{}_Ntau{}_dt{}.inp'.format(ns,U,Nt,Ntau,dt)
    else:
        flin = input_file
    GenInputFile(flin,sysparams,solverparams)
    prog = runpath + '/exe/gw.x'

    log_flag = ''
    if len(log_file) > 0:
        log_flag = ' > ' + log_file
        
    os.system("mpirun " + prog + ' ' + flin + ' ' + output_file + log_flag)
#----------------------------------------------------------------------     

if __name__ == '__main__':

    os.chdir("build")
    os.system("../osx_cmake.sh")
    os.chdir("../")

    Nk=10
    thop=1.0
    U=2.0
    V=0.0
    mu=1.0
    beta=10.0
    sysparams = GenSysParams(Nk,thop,U,V,mu,beta)

    Nt = 200
    Ntau = 700
    dt = 0.01
    solverparams = GenSolverParams(Nt,Ntau,dt,MatsMaxErr= 1.0e-10,BootstrapMaxErr=1.0e-10)

    #fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)

    color_palette = ['blue', 'orange', 'green', 'red']

    #ax2.set_xlabel(r'$t$')
    #ax2.set_ylabel('energy')
    #ax1.set_ylabel(r'$n_1(t)$')
    
    output_file = 'out/gw_'
    Run(sysparams,solverparams,output_file,runpath='./')

    data = np.loadtxt(output_file + '_occupation.dat')
    t = data[:,0]        
    #ax1.plot(t,data[:,1],label=approx,c=color_palette[i])

    t,Ekin,Epot,Etot = np.loadtxt(output_file + '_energy.dat',unpack=True)
    #ax2.plot(t,Ekin,c=color_palette[i])
    #ax2.plot(t,Etot,c=color_palette[i],linestyle='--')

    #ax1.legend(loc='upper right',frameon=False)
    #plt.show()

