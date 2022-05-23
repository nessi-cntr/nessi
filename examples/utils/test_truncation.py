import sys
import os
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
sys.path.insert(0, '../libcntr/python/')
sys.path.insert(0, '../libcntr/python3/')
from ReadCNTR import write_input_file
#----------------------------------------------------------------------
def GenInputFile(input_file,u0,u1,neps,max_energy,min_energy,nt,beta,
                 ntau,h,tmax,tc):
    inp = {
        'u0':u0,
        'u1':u1,
        'neps':neps,
        'max_energy':max_energy,
        'min_energy':min_energy,
        'nt': nt,
        'beta':beta,
        'ntau':ntau,
        'h':h,
        'tmax':tmax,
        'tc':tc
        }
    write_input_file(inp, input_file)
#----------------------------------------------------------------------
def RunTestTruncation(u0,u1,neps,max_energy,min_energy,nt,beta,
                       ntau,h,tmax,tc,output_file,runpath='./'):  
    flin = './inp/param_test_truncation.inp'
    GenInputFile(flin,u0,u1,neps,max_energy,min_energy,nt,beta,
                 ntau,h,tmax,tc)
    
    prog = runpath + 'exe/test_truncation_bethe_quench.x'
    os.system(prog + ' ' + flin + ' ' + output_file)
#----------------------------------------------------------------------     
if __name__ == '__main__':
    input_file = './inp/param_test_truncation.inp'
    output_file = './out/'

    u0=0
    u1=1
    neps=20
    max_energy=0.4
    min_energy=0.0
    nt=500
    beta=100
    ntau = 1000
    h=0.04
    tmax = 10000
    tc=500#1000
    
  
    #RunTestTruncation(u0,u1,neps,max_energy,min_energy,nt,beta,
    #                   ntau,h,tmax,tc,output_file,runpath='./')

#---Plotting-----------------------------------------------------------
    cmap=matplotlib.cm.get_cmap('jet')
    C_LIST = cmap(np.linspace(0, 1, neps))
    fig1, ax1 = plt.subplots(1, 1)
    color_ind=0
    for j in range(neps):
        filename=output_file+"temp_nk_tc_"+str(tc)+"_h_"+str(h)+"_ui_"+str(u0)+".0_uf_"+str(u1)+".00_k_"+str(j)+"_tmax_"+str(tmax)+".out"
        ek=np.loadtxt(filename,usecols=(3,))
        t=np.loadtxt(filename,usecols=(5,))
        nk=np.loadtxt(filename,usecols=(7,))
        nk_trunc=np.loadtxt(filename,usecols=(9,))
        line="-"
        color=C_LIST[color_ind]
        label=r'$\epsilon=$'+str(ek[0])
        if (j<=1 or j==neps-1):
            ax1.semilogx(t*h,nk_trunc,color=color,ls=line,label=label)
            ax1.semilogx(tc*h,nk_trunc[tc],color=color,ls=None,marker='x')
        else:
            ax1.semilogx(t*h,nk_trunc,color=color,ls=line)
            ax1.semilogx(tc*h,nk_trunc[tc],color=color,ls=None,marker='x')
        color_ind+=1
    ##ADD CBAR FOR ENERGY INDEX??
    ax1.set_xlim(h,tmax*h)
    ax1.set_ylim(0,0.59)
    ax1.set_ylabel(r"$n_{\epsilon}(t)$")
    ax1.set_xlabel(r"$t$")

    fig1.subplots_adjust(top=0.965,
                         bottom=0.11,
                         left=0.1,
                         right=0.965,
                         hspace=0.0,
                         wspace=0.2)


    fig1.savefig(output_file+"Truncation_test.png")
    plt.show()
