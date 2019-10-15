import sys
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from ReadCNTR import write_input_file
#----------------------------------------------------------------------
def GenInputFile(input_file,Ntau,SolveOrder):
    inp = {'Ntau': Ntau, 'SolveOrder': SolveOrder}
    write_input_file(inp, input_file)
#----------------------------------------------------------------------
def RunTestEquilibrium(Ntau,SolveOrder,output_file,runpath='./'):  
    flin = './inp/test_equilibrium.inp'
    GenInputFile(flin,Ntau,SolveOrder)
    
    prog = runpath + 'exe/test_equilibrium.x'
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
    ptau = np.linspace(1.0,3.0,20);
    Ntau_vals = np.int_(np.ceil(np.power(10,ptau)));

    output_file = 'out/test_equilibrium.dat'
    
    if os.path.exists(output_file):
        os.remove((output_file))
              
    fig,ax = plt.subplots()

    for Ntau in Ntau_vals:
        Ntau_even = Ntau
        if not Ntau % 2 == 0:
            Ntau_even = Ntau + 1           
        RunTestEquilibrium(Ntau_even,SolveOrder,output_file)

    err1,err2 = np.loadtxt(output_file,unpack=True)
    logn = np.log10(Ntau_vals)
    log_err1 = np.log10(err1)
    log_err2 = np.log10(err2)

    ax.scatter(logn,log_err1,label='Fourier',marker='o',s=100,edgecolor='k',facecolor='lightgray')
    ax.scatter(logn,log_err2,label='fix point',marker='s',s=100,edgecolor='r',facecolor='mistyrose')
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(logn[6:],log_err1[6:])
    print('order(Fourier) : ' + '%10.2f' % slope)
    ax.plot(logn,slope*logn+intercept,'k--',label=r'$\sim \mathcal{O}(h_\tau^{' + '%10.2f' % -slope + '})$')
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(logn[6:],log_err2[6:])
    print('order(fix point) : ' + '%10.2f' % slope)
    ax.plot(logn,slope*logn+intercept,'r--',label=r'$\sim \mathcal{O}(h_\tau^{' + '%10.2f' % -slope + '})$' )
    

    ax.set_title(r'$k={}$'.format(SolveOrder),fontsize=20)
    ax.legend(loc='lower left',frameon=False,fontsize=16)
    ax.set_xlabel(r'$\log(N_\tau)$')
    ax.set_ylabel(r'$\log(\mathrm{err})$')
    plt.show()

