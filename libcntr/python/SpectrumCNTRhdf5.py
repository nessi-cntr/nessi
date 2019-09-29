# ----------------------------------------------------------------------

import h5py
import numpy as np
import scipy as sp
from scipy import interpolate,integrate

import sys
import matplotlib.pyplot as plt

from operator import itemgetter
from scipy.integrate import simps

# ----------------------------------------------------------------------
# Read G^R(tstp,t) from 'G_slice' for t<tstp
# Outputs are given and G^R(t') = G^R(tstp,tstp-t') and list of t'
# [note]; t_rel is given every dt
def read_gf_ret_slice(filename,name_green,dt,tstp=0):
    
    fd = h5py.File(filename)
    name_tstp = 't{}'.format(tstp)
    G_ret_ = fd[name_green][name_tstp]['ret'][:,:,:]
    Nt = len(G_ret_)-1
    t_rel = np.linspace(0.0,float(Nt)*dt,Nt+1)
    G_ret = np.flip(G_ret_,0)

    return t_rel, G_ret

# ----------------------------------------------------------------------
# Read G^R(t_av;t_rel) and G^<(t_av;t_rel) from 'G_tavtrel' for specified tav and trel>=0
# [note]; t_rel is given every 2dt
def read_gf_ret_let_tavtrel(filename,name_green,dt,tav=0):

    fd = h5py.File(filename)
    name_tav = '{}'.format(tav)
    G_ret = fd[name_green]['gtr'][name_tav][:,:,:]-fd[name_green]['les'][name_tav][:,:,:]
    G_les = fd[name_green]['les'][name_tav][:,:,:]

    Nt = len(G_ret)-1
    t_rel = np.linspace(0.0,float(Nt)*2.0*dt,Nt+1)
    
    return t_rel, G_ret, G_les

# ----------------------------------------------------------------------
def complex_quad(f_t,t0,t1):
    
    xj = 1.0j
    
    def re_f_t(t):
        return np.real(f_t(t))
    
    def im_f_t(t):
        return np.imag(f_t(t))
    
    Re_F = integrate.quad(re_f_t,t0,t1)
    Im_F = integrate.quad(im_f_t,t0,t1)
    return Re_F[0]+xj*Im_F[0], Re_F[1]+xj*Im_F[1]

# ----------------------------------------------------------------------
#  f(w)=\int^tmax_t_min dt v(t)e^{i\omega t}
def fourier_t2w(v_t,tgrd,wgrd,method):
    
    methods = ['linear','quadratic','cubic','simpson']
    assert (method in methods), '@fourier_t2w. This method does not exist'
    
    xj = 1.0j
    Nt = len(v_t)-1
    dt = tgrd[1]-tgrd[0]
    
    v_w = np.zeros(len(wgrd),dtype='complex')

    if method == 'linear':
        v_w = linear_fourier_t2w(v_t,tgrd,wgrd)

    elif method == 'quadratic' or method == 'cubic':
 
        f_t = interpolate.interp1d(tgrd,v_t, kind=method, bounds_error=False,fill_value=0.0)

        for iw in range(len(wgrd)):
            def f_t_exp_iwt(t):
                return f_t(t)*np.exp(xj*wgrd[iw]*t)

            v_w[iw],error = complex_quad(f_t_exp_iwt,tgrd[0],tgrd[-1])
    
    elif method == 'simpson':
        for iw in range(len(wgrd)):
            f_re = np.cos(wgrd[iw]*tgrd[:])*v_t[:]
            f_im = np.sin(wgrd[iw]*tgrd[:])*v_t[:]
            v_w_re = simps(f_re,dx=dt)
            v_w_im = simps(f_im,dx=dt)
            v_w[iw] = v_w_re + xj*v_w_im
                
    return v_w

#----------------------------------------------------------------------
def linear_fourier_t2w(v_t,tgrd,wgrd):

    xj = 1.0j

    Nt = len(v_t)-1
    dt = tgrd[1]-tgrd[0]
    v_w = np.zeros(len(wgrd),dtype='complex')
    
    for i in range(len(wgrd)):
        ome_ = wgrd[i]
    
        if ome_ == 0.0:
            v_w[i] += 0.5*v_t[0]*dt
            for it in range(1,Nt_):
                v_w += v_t[it]*dt
            v_w[i] += 0.5*v_t[Nt]*dt
        else:
            for it in range(0,Nt):
                v_w[i] += v_t[it]*np.exp(xj*ome_*tgrd[it])
            
            v_w[i] *= 4.0*np.sin(ome_*dt/2.0)*np.sin(ome_*dt/2.0)/ome_/ome_/dt
            
            v_w[i] += (v_t[Nt]*np.exp(xj*ome_*tgrd[-1])-v_t[0]*np.exp(xj*ome_*tgrd[0]))\
                *(1.0/xj/ome_-(np.exp(-xj*ome_*dt)-1.0)/ome_/ome_/dt)
    return v_w
#----------------------------------------------------------------------
def evaluate_spectrum(green_ret,tgrd,wgrd,method):
    
    w_size = len(wgrd)
    element_size = green_ret.shape[1]

    Aret = np.zeros(shape=(w_size,element_size))

    for ie in range(element_size):
        green_ret_damp_ = green_ret[:,ie,ie]
        green_ret_w_ = fourier_t2w(green_ret_damp_,tgrd,wgrd,method)
        Aret[:,ie] = -1.0/np.pi*np.imag(green_ret_w_[:])

    return Aret

#----------------------------------------------------------------------
def evaluate_spectrum_windowed(green_ret,tgrd,wgrd,method,Fwindow):
    
    w_size = len(wgrd)
    element_size = green_ret.shape[1]
    
    Aret = np.zeros(shape=(w_size,element_size))
    
    for ie in range(element_size):
        green_ret_damp_ = green_ret[:,ie,ie]*Fwindow(tgrd[:])
        green_ret_w_ = fourier_t2w(green_ret_damp_,tgrd,wgrd,method)
        Aret[:,ie] = -1.0/np.pi*np.imag(green_ret_w_[:])
    
    return Aret

#----------------------------------------------------------------------
def evaluate_occupation(green_less,tgrd,wgrd,method):
    
    w_size = len(wgrd)
    element_size = green_less.shape[1]
    
    Nles = np.zeros(shape=(w_size,element_size))
    
    for ie in range(element_size):
        green_less_damp_ = green_less[:,ie,ie]
        green_less_w_ = fourier_t2w(green_less_damp_,tgrd,wgrd,method)
        if tgrd[0] == 0.0: # for the case only t>=t' for G^< is given
            Nles[:,ie] = 1.0/(2.0*np.pi)*np.imag(green_less_w_[:]-np.conj(green_less_w_[:]))
        else:
            Nles[:,ie] = 1.0/(2.0*np.pi)*np.imag(green_less_w_[:])
    
    return Nles
#----------------------------------------------------------------------
def evaluate_occupation_windowed(green_less,tgrd,wgrd,method,Fwindow):
    
    w_size = len(wgrd)
    element_size = green_less.shape[1]
    
    Nles = np.zeros(shape=(w_size,element_size))
    
    for ie in range(element_size):
        green_less_damp_ = green_less[:,ie,ie]*Fwindow(tgrd[:])
        green_less_w_ = fourier_t2w(green_less_damp_,tgrd,wgrd,method)
        if tgrd[0] == 0.0: # for the case only t>=t' for G^< is given
            Nles[:,ie] = 1.0/(2.0*np.pi)*np.imag(green_less_w_[:]-np.conj(green_less_w_[:]))
        else:
            Nles[:,ie] = 1.0/(2.0*np.pi)*np.imag(green_less_w_[:])
    
    return Nles
# ----------------------------------------------------------------------
'''
if __name__ == "__main__":

    filename = sys.argv[1]
    dt = 0.03
    name_green = 'Gloc_slice'
    t_rel1,G_ret1 = read_gf_ret_slice(filename,name_green,dt,2000)
    name_green = 'Gloc_tavrel'
    t_rel,G_ret,G_les = read_gf_ret_let_tavtrel(filename,name_green,dt,1000)
    
    
    ome = np.linspace(-4,4,800)
    method = 'linear'
    def Fwindow(t):
        return np.exp(-t*t/600.0)
    
    A_w = evaluate_spectrum_windowed(G_ret,t_rel,ome,method,Fwindow)
    N_w = evaluate_occupation_windowed(G_les,t_rel,ome,method,Fwindow)
    
    #v_w1 = fourier_t2w(G_ret1[:,0,0],t_rel1-10.0,ome,'linear')
    #v_w3 = fourier_t2w(G_ret1[:,0,0],t_rel1-10.0,ome,'cubic')
    #v_w2 = fourier_t2w(G_ret1[:,0,0],t_rel1-10.0,ome,'simpson')
    
    fig,ax = plt.subplots()
    #ax.plot(ome,np.imag(v_w1))
    #ax.plot(ome,np.imag(v_w2))
    #ax.plot(ome,np.imag(v_w3))
    ax.plot(ome,A_w[:,0])
    ax.plot(ome,A_w[:,0]/(np.exp(40.0*ome[:])+1))
    ax.plot(ome,N_w[:,1])
    
    plt.show()
'''
