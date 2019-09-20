import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
# ---- NESSY python modules ----
from ReadCNTR import write_input_file
from SpectrumCNTRhdf5 import read_gf_ret_let_tavtrel,evaluate_spectrum_windowed
#----------------------------------------------------------------------
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z
#----------------------------------------------------------------------
def GenSysParams(Hopping,El_Ph_g,Phfreq_w0,MuChem_MF,Beta):
    sysparams = {'Hopping': Hopping,
                 'El_Ph_g': El_Ph_g,
                 'Phfreq_w0': Phfreq_w0,
                 'MuChem_MF': MuChem_MF,
                 'beta': Beta
        }
    return sysparams

#----------------------------------------------------------------------
def Sin2_env(time,tend):
    dvalue = 0.0
    if time>0.0 and time < tend:
        x =time*np.pi/tend
        dvalue = np.sin(x)*np.sin(x)
    return dvalue

#----------------------------------------------------------------------
def GenExciteParams(Excite_ElPh_g,Excite_Hop,Excite_Freq,Tend,Nt,dt,runpath='./'):
    
    flin_dHopping = runpath+'inp/dHopping_field.txt'
    flin_dElPh_g = runpath+'inp/dElPh_g_field.txt'
    
    ##
    exciteparams = {'dHopping': '--' + flin_dHopping,
        'dEl_Ph_g': '--' + flin_dElPh_g 
    }
    
    time = np.linspace(-dt,float(Nt)*dt,Nt+2)
    dElPh_g = np.array([])
    dHopping = np.array([])
    dElPh_g  = np.append(dElPh_g,0.0)
    dHopping = np.append(dHopping,0.0)
    
    for it in range(Nt+1):
        dElPh_g_ = Excite_ElPh_g*Sin2_env(time[it+1],Tend)*np.cos(time[it+1]*Excite_Freq)
        dElPh_g = np.append(dElPh_g,dElPh_g_)
        dHopping_ = Excite_Hop*Sin2_env(time[it+1],Tend)*np.cos(time[it+1]*Excite_Freq)
        dHopping = np.append(dHopping,dHopping_)
    
    np.savetxt(flin_dElPh_g,dElPh_g,fmt='%.10f')
    np.savetxt(flin_dHopping,dHopping,fmt='%.10f')
    
    return exciteparams
#----------------------------------------------------------------------
def GenSolverParams(Nt,Ntau,dt,SaveGreen=False,MatsMaxIter=800,
                        MatsMaxErr= 1.0e-7,BootstrapMaxIter= 20,
                        BootstrapMaxErr=1.0e-7,CorrectorSteps=5,SolverOrder=5):

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
                    'CorrectorSteps': CorrectorSteps,
                    'SolverOrder': SolverOrder,
                    'SaveGreen': savegf
                        }
    return solverparams

def GenOutParams(OutEvery):
    
    outparams = { 'OutEvery':OutEvery }

    return outparams

#----------------------------------------------------------------------
def GenInputFile(input_file,sysparams,exciteparams,solverparams,outparams):
    #inp = {**sysparams, **exciteparams, **solverparams, **outparams}
    inp = merge_two_dicts(sysparams, exciteparams)
    inp = merge_two_dicts(inp, solverparams)
    inp = merge_two_dicts(inp, outparams)
    write_input_file(inp, input_file)
#----------------------------------------------------------------------
def RunHolstein(sysparams,exciteparams,solverparams,outparams,approx_sigma,output_file,
                    log_file='',runpath='./'):

    Nt = solverparams['Nt']
    Ntau = solverparams['Ntau']
    dt = solverparams['dt']
    flin = runpath+'/inp/Holstein_Nt{}_Ntau{}_dt{}.inp'.format(Nt,Ntau,dt)
    
    GenInputFile(flin,sysparams,exciteparams,solverparams,outparams)
    if approx_sigma == 'sMig':
        prog = runpath + 'exe/Holstein_bethe_Migdal.x'
    elif approx_sigma == 'uMig':
        prog = runpath + 'exe/Holstein_bethe_uMig.x'
    else:
        print('[ERROR] Self-energy approximation not recognized.')
        exit()

    log_flag = ''
    if len(log_file) > 0:
        log_flag = ' > ' + log_file

    os.system(prog + ' ' + flin + ' ' + output_file + log_flag)
#----------------------------------------------------------------------     
if __name__ == '__main__':
    
    # Parameter of the Holstein model
    thop=1.0
    g=0.5
    w0=1.0
    mu_mf=0.5
    beta=10.0
    sysparams = GenSysParams(thop,g,w0,mu_mf,beta)

    # calculation parameters 
    Nt = 400
    Ntau = 500
    dt = 0.03
    solverparams = GenSolverParams(Nt,Ntau,dt)
    
    # excitation parameters
    Excite_ElPh_g = 0.0
    Excite_Hop = 1.0
    Excite_Freq = 0.0
    Excite_End = 0.6
    exciteparams = GenExciteParams(Excite_ElPh_g,Excite_Hop,Excite_Freq,Excite_End,Nt,dt)

    #  output parameters
    OutEvery = int(Nt/10)
    outparams = GenOutParams(OutEvery)
    tstp_plot = int(Nt/2) #should be integer multiplication of OutEvery
    

    # impurity solver
    sigma_approx = ['uMig']
    
    #plot 
    fig,ax = plt.subplots(3,1,sharex=True)
    fig2,ax2 = plt.subplots(1,2)
    ax[-1].set_xlabel(r'$t$')
    ax[0].set_ylabel(r'$n(t)$')
    ax[1].set_ylabel(r'$X(t),P(t)$')
    ax[2].set_ylabel(r'Energy')
    ax2[0].set_xlabel(r'$\omega$')
    ax2[1].set_xlabel(r'$\omega$')
    ax2[0].set_ylabel(r'$A(\omega)$')
    ax2[1].set_ylabel(r'$B(\omega)$')
        
    for approx in sigma_approx:
        output_file = './out/Holstein_' + approx
        RunHolstein(sysparams,exciteparams,solverparams,outparams,approx,output_file)
        
        #Plotting observables
        data = np.loadtxt(output_file + '_obs.dat')
        t = data[:,0]
        ax[0].plot(t,2.0*data[:,1],label=approx+r':$n_{\rm tot}$')
        ax[1].plot(t,data[:,3],label=approx+r':$X$')
        ax[1].plot(t,data[:,4],label=approx+r':$P$')
        ax[2].plot(t,data[:,6],label=approx+r'$:E_{\rm kin}$')
        ax[2].plot(t,data[:,7]+data[:,8],label=approx+r'$:E_{\rm nx}$')
        ax[2].plot(t,data[:,9]+data[:,10],label=approx+r'$:E_{\rm ph}$')
        ax[2].plot(t,data[:,11],label=approx+r'$:E_{\rm tot}$')

        #Plotting Spectral function
        ome = np.linspace(-4,4,800)
        method = 'linear'
        def Fwindow(t):
            return np.exp(-t*t/400.0)
        
        t_rel,G_ret,G_les = read_gf_ret_let_tavtrel(output_file + '_green.h5','Gloc_tavrel',dt,tstp_plot)
        t_rel,D_ret,D_les = read_gf_ret_let_tavtrel(output_file + '_green.h5','Dloc_tavrel',dt,tstp_plot)

        A_w = evaluate_spectrum_windowed(G_ret,t_rel,ome,method,Fwindow)
        B_w = evaluate_spectrum_windowed(D_ret,t_rel,ome,method,Fwindow)
        
        ax2[0].plot(ome,A_w[:,0])
        ax2[1].plot(ome,B_w[:,0])
        
    ax[0].legend(loc='upper right',frameon=False)
    ax[1].legend(loc='upper right',frameon=False)
    ax[2].legend(loc='upper right',frameon=False)

    plt.subplots_adjust(wspace=0.55)
    
    ax2[0].legend(loc='upper right',frameon=False)
    ax2[1].legend(loc='upper right',frameon=False)
    plt.show()

