import sys
import os
import numpy as np
from scipy.integrate import simps
from scipy import stats
import matplotlib.pyplot as plt
from ReadCNTR import write_input_file
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def Iexact(x):
    return -1j*(np.exp(1j*x)-1.0)
#----------------------------------------------------------------------
def RunIntegration(npts,h,k,output_file,input_file='',log_file='',runpath='./'):
    prog = runpath + '/exe/integration.x'
    inp = {'npts': npts, 'h': h, 'k': k}

    if len(input_file) == 0:
        flin = './inp/integ.inp'
    else:
        flin = input_file
        
    write_input_file(inp, flin)
        
    log_flag = ''
    if len(log_file) > 0:
        log_flag = ' > ' + log_file
        
    os.system(prog + ' ' + flin + ' ' + output_file + log_flag)
#----------------------------------------------------------------------     

if __name__ == '__main__':

    xmax = 2.5*np.pi
    npts_vals = [2**j for j in range(4,10)]
    logn = np.log10(npts_vals)

    err_greg = np.zeros((len(npts_vals),6))
    err_simps = np.zeros(len(npts_vals))

    for ir,npts in enumerate(npts_vals):
        h = xmax/npts
        for k in range(1,6):
            flout = 'out/integ_k{}.dat'.format(k)
            RunIntegration(npts,h,k,flout,runpath='./')

            x,err = np.loadtxt(flout,unpack=True)
            err_greg[ir,k-1] = np.mean(err)
        for n in range(1,npts+1):
            xr = h*np.arange(0,n+1)
            c_int = simps(np.cos(xr),dx=h)
            s_int = simps(np.sin(xr),dx=h)
            err_simps[ir] += np.abs(c_int + 1.j*s_int -Iexact(xr[n]))/npts

    fig,ax = plt.subplots(figsize=(5,4))

    ax.set_xlabel(r'log($N$)')
    ax.set_ylabel('log(err.)')

    palette = ['blue','green','orange','red','purple']
    
    for k in range(1,6):
        slope, intercept, r_value, p_value, std_err = stats.linregress(logn,np.log10(err_greg[:,k-1]))
        print('expected order: {} | approx. order: {:10.2f}'.format(k+2,-slope))
        ax.plot(logn,slope*logn+intercept,c=palette[k-1],ls='--',label=r'$p={:10.2f}$'.format(-slope))
        ax.scatter(logn,np.log10(err_greg[:,k-1]),edgecolor='k',facecolor = palette[k-1],s=80)

    slope, intercept, r_value, p_value, std_err = stats.linregress(logn,np.log10(err_simps))
    ax.plot(logn,slope*logn+intercept,c='gray',ls='--',label='simps.')
    ax.scatter(logn,np.log10(err_simps),edgecolor='k',facecolor ='gray',s=80)


    ax.legend(loc='best',ncol=2,fontsize=12,frameon=False)

    plt.show()

    
