import sys
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from ReadCNTR import write_input_file
#----------------------------------------------------------------------
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z
#----------------------------------------------------------------------
def GenField(E0,omega,Np,Nt,h,file_field):
    ts = np.linspace(0.0,Nt*h,Nt+1)
    t0 = 2.0*np.pi/(omega*Np)
    Ef = E0 * np.exp(-4.6*(ts-t0)**2/t0**2)*np.sin(omega*(ts-t0))
    Ef = np.insert(Ef, 0, 0.0)
    np.savetxt(file_field,Ef,fmt='%.10f')
#----------------------------------------------------------------------
def GenSysParams(Nk,HoppingT,HubbardU,V,MuChem,Beta):
    sysparams = {'Nk': Nk,
                 'HoppingT': HoppingT,
                 'HubbardU': HubbardU,
                 'V': V,
                 'MuChem': MuChem,
                 'beta': Beta
        }
    return sysparams
#----------------------------------------------------------------------
def GenSolverParams(Nt,Ntau,h,SaveGreen=False,SaveMomentum=False,MatsMaxIter=100,
                        MatsMaxErr= 1.0e-8,
                        BootstrapMaxIter= 20,BootstrapMaxErr=1.0e-8,
                        TimeMaxErr=1.0e-8,CorrectorSteps=10,output=100):

    savegf = 0
    if SaveGreen:
        savegf = 1

    savegk = 0
    if SaveMomentum:
        savegk = 1
    
    solverparams = {'Nt': Nt,
                    'Ntau': Ntau,
                    'h': h,
                    'MatsMaxIter': MatsMaxIter,
                    'MatsMaxErr': MatsMaxErr,
                    'BootstrapMaxIter': BootstrapMaxIter,
                    'BootstrapMaxErr': BootstrapMaxErr,
                    'TimeMaxErr': TimeMaxErr,
                    'CorrectorSteps': CorrectorSteps,
                    'SaveGreen': savegf,
                    'SaveMomentum': savegk
                        }
    if SaveMomentum:
        solverparams.update({'output':output})

    return solverparams
#----------------------------------------------------------------------
def GenInputFile(input_file,sysparams,solverparams,file_field):
    inp = sysparams
    inp = merge_two_dicts(inp, solverparams)
    inp = merge_two_dicts(inp, {'Epulse': '--'+file_field})
    write_input_file(inp, input_file)
#----------------------------------------------------------------------
def Run(sysparams,solverparams,file_field,output_file,mpicmd,
                        input_file='',log_file='',runpath=''):
    Nk = sysparams['Nk']
    U = sysparams['HubbardU']
    V = sysparams['V']
    Nt = solverparams['Nt']
    Ntau = solverparams['Ntau']
    h = solverparams['h']
    if len(input_file) == 0:
        flin = './inp/gw_Nk{}_U{}_V{}_Nt{}_Ntau{}_h{}.inp'.format(Nk,U,V,Nt,Ntau,h)
    else:
        flin = input_file
    GenInputFile(flin,sysparams,solverparams,file_field)
    prog = runpath + '/exe/gw.x'

    log_flag = ''
    if len(log_file) > 0:
        log_flag = ' > ' + log_file
        
    os.system(mpicmd + ' ' + prog + ' ' + flin + ' ' + log_flag)
#----------------------------------------------------------------------     
def ReadData(fname):
	f = h5py.File(fname, 'r')

	Nt = f['parm']['Nt'][0]
	dt = f['parm']['dt'][0]
	Epulse = np.array(np.real(f['parm']['Epulse']['data'][1:,0,0]))
	Ekin = np.array(np.real(f['obs']['Ekin']['data'][1:,0,0]))
	Epot = np.array(np.real(f['obs']['Epot']['data'][1:,0,0]))

	f.close()

	ts = np.linspace(0.0,Nt*dt,Nt+1)

	return ts, Epulse,Ekin,Epot
#---------------------------------------------------------------------- 
if __name__ == '__main__':

    Nk=10
    thop=1.0
    U=2.0
    V=0.0
    mu=1.0
    beta=10.0
    sysparams = GenSysParams(Nk,thop,U,V,mu,beta)

    Nt = 1
    Ntau = 100
    h = 0.02
    solverparams = GenSolverParams(Nt,Ntau,h,SaveGreen=True,SaveMomentum=True,MatsMaxErr= 1.0e-10,BootstrapMaxErr=1.0e-10,output=1)

    # field 
    E0 = 0.2
    Np = 1
    omega = 1.0
    file_field = 'inp/Epulse_E{}_Np{}_w{}.txt'.format(E0,Np,omega)
    GenField(E0,omega,Np,Nt,h,file_field)
    
    output_file = 'out/gw_'
    mpicmd = 'mpirun -n 2'
    Run(sysparams,solverparams,file_field,output_file,mpicmd,runpath='./')

    # ---- read from file ----
    fname = 'out/data_gw.h5'
    ts, Epulse, Ekin, Epot = ReadData(fname)

    fig, ax = plt.subplots(2,1,sharex=True)
    ax[0].plot(ts, Epulse, c='blue')
    ax[1].plot(ts, Ekin - Ekin[0], c='k')

    ax[0].set_ylabel(r'$E(t)$')
    ax[1].set_ylabel(r'$E_\mathrm{kin}(t) - E_\mathrm{kin}(0) $')
    ax[1].set_xlabel('time')

    plt.show()

