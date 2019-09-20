import sys
import os
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
from ReadCNTR import write_input_file
#----------------------------------------------------------------------

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

    npts = 50
    h = 0.05*np.pi

    fig,ax = plt.subplots()
    ax.set_yscale("log", nonposy='clip')        

    ax.set_xlabel(r'$x/\pi$')
    ax.set_ylabel('abs. error')

    palette = ['blue','green','orange','red','purple']
    
    for k in range(1,6):
        flout = 'out/integ_k{}.dat'.format(k)
        RunIntegration(npts,h,k,flout,runpath='./')

        x,err = np.loadtxt(flout,unpack=True)

        ax.plot(x[1:]/(np.pi),err[1:],c = palette[k-1],label=r'$k={}$'.format(k),marker='o',markersize=4)

    err_simps = np.zeros(npts+1)
    for n in range(1,npts+1):
        xr = h*np.arange(0,n+1)
        c_int = simps(np.cos(xr),dx=h)
        s_int = simps(np.sin(xr),dx=h)
        err_simps[n] = np.abs(c_int + 1.j*s_int + 1.j*(np.exp(1.j*xr[n])-1.0)) 
 
    ax.plot(x[1:]/(np.pi),err_simps[1:],c='gray',label='simps.')
    ax.set_ylim(1.0e-10,1.0e-3)

    ax.legend(loc='upper right',ncol=3,fontsize=12)
    plt.show()

    
