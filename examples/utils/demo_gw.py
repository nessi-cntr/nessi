import sys
import os
import glob
import numpy as np
from ReadCNTR import write_input_file
#----------------------------------------------------------------------
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z
#----------------------------------------------------------------------
def GenField(E0,omega,Np,Nt,dt,file_field):
    ts = np.linspace(0.0,Nt*dt,Nt+1)
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
def GenSolverParams(Nt,Ntau,dt,SaveGreen=False,SaveMomentum=False,MatsMaxIter=100,
                        MatsMaxErr= 1.0e-8,
                        BootstrapMaxIter= 20,BootstrapMaxErr=1.0e-8,
                        TimeMaxErr=1.0e-8,CorrectorSteps=10):

    savegf = 0
    if SaveGreen:
        savegf = 1

    savegk = 0
    if SaveMomentum:
        savegk = 1
    
    solverparams = {'Nt': Nt,
                    'Ntau': Ntau,
                    'dt': dt,
                    'MatsMaxIter': MatsMaxIter,
                    'MatsMaxErr': MatsMaxErr,
                    'BootstrapMaxIter': BootstrapMaxIter,
                    'BootstrapMaxErr': BootstrapMaxErr,
                    'TimeMaxErr': TimeMaxErr,
                    'CorrectorSteps': CorrectorSteps,
                    'SaveGreen': savegf,
                    'SaveMomentum': savegk
                        }
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
    dt = solverparams['dt']
    if len(input_file) == 0:
        flin = './inp/gw_Nk{}_U{}_V{}_Nt{}_Ntau{}_dt{}.inp'.format(Nk,U,V,Nt,Ntau,dt)
    else:
        flin = input_file
    GenInputFile(flin,sysparams,solverparams,file_field)
    prog = runpath + '/exe/gw.x'

    log_flag = ''
    if len(log_file) > 0:
        log_flag = ' > ' + log_file
        
    os.system(mpicmd + ' ' + prog + ' ' + flin + ' ' + output_file + log_flag)
#----------------------------------------------------------------------     

if __name__ == '__main__':

    Nk=10
    thop=1.0
    U=2.0
    V=0.0
    mu=1.0
    beta=10.0
    sysparams = GenSysParams(Nk,thop,U,V,mu,beta)

    Nt = 200
    Ntau = 200
    dt = 0.01
    solverparams = GenSolverParams(Nt,Ntau,dt,MatsMaxErr= 1.0e-10,BootstrapMaxErr=1.0e-10)

    # field 
    E0 = 0.0
    Np = 1
    omega = 1.0
    file_field = 'inp/Epulse_E{}_Np{}_w{}.txt'.format(E0,Np,omega)
    GenField(E0,omega,Np,Nt,dt,file_field)

    color_palette = ['blue', 'orange', 'green', 'red']
    
    output_file = 'out/gw_'
    mpicmd = 'mpirun -n 2'
    Run(sysparams,solverparams,file_field,output_file,mpicmd,runpath='./')

