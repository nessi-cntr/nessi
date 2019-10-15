import sys
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from ReadCNTR import write_input_file
#----------------------------------------------------------------------
def GenInputFile(input_file,Nt,Tmax,Ntau,SolveOrder):
    inp = {
        'Nt': Nt,
        'Tmax': Tmax,        
        'Ntau': Ntau,
        'SolveOrder': SolveOrder
        }
    write_input_file(inp, input_file)
#----------------------------------------------------------------------
def RunTestEquilibrium(Nt,Tmax,Ntau,SolveOrder,output_file,runpath='./'):  
    flin = './inp/test_nonequilibrium.inp'
    GenInputFile(flin,Nt,Tmax,Ntau,SolveOrder)
    
    prog = runpath + 'exe/test_nonequilibrium.x'
    os.system(prog + ' ' + flin + ' ' + output_file)
#----------------------------------------------------------------------     

if __name__ == '__main__':
    if len(sys.argv) > 1:
        SolveOrder = int(sys.argv[1])
    else:
        print("As first argument provide the solve order k=1,...,5")
        exit()
    if SolveOrder > 5 :
        print("As first argument provide the solve order k=1,...,5")
        exit()
    input_file = './inp/test_nonequilibrium.inp'
    output_file = './out/test_nonequilibrium.dat'
    if os.path.isfile(output_file):
        os.remove(output_file)
       
    Ntau = 800
    Tmax = 5.0

    nt_val = np.power(10,np.arange(start=1.0,stop=2.0,step=0.05))

    fig,ax = plt.subplots()
    
    for Nt in nt_val:
        RunTestEquilibrium(Nt,Tmax,Ntau,SolveOrder,output_file,runpath='./')
        
    err_dyson,err_vie2 = np.loadtxt(output_file,unpack=True)
    log_nt = np.log10(nt_val)
    log_err_dyson = np.log10(err_dyson)
    log_err_vie2 = np.log10(err_vie2)

    ax.scatter(log_nt,log_err_dyson,marker='o',s=100,edgecolor='k',facecolor='lightgray',\
                   label='dyson (k = {})'.format(SolveOrder))          
    ax.scatter(log_nt,log_err_vie2,marker='s',s=100,edgecolor='r',facecolor='mistyrose',\
                   label='vie2 (k = {})'.format(SolveOrder))

    slope, intercept, r_value, p_value, std_err = stats.linregress(log_nt[10:],log_err_dyson[10:])
    print('order(Dyson) : ' + '%10.2f' % slope)
    ax.plot(log_nt,slope*log_nt+intercept,'k--',label=r'$\sim \mathcal{O}(h^{' + '%10.2f' % -slope + '})$')
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_nt[10:],log_err_vie2[10:])
    print('order(VIE2) : ' + '%10.2f' % slope)
    ax.plot(log_nt,slope*log_nt+intercept,'r--',label=r'$\sim \mathcal{O}(h^{' + '%10.2f' % -slope + '})$' )
    
    ax.legend(loc='lower left',frameon=False,fontsize=16)
    ax.set_xlabel(r'$\log(N_t)$')
    ax.set_ylabel(r'$\log(\mathrm{err})$')
    plt.show()

