#!/usr/bin/env python
# coding: utf-8

# # Demo: Steady state solution of the Anderson impurity model

# ## Overview
# 
# This notebook provides a perturbative solution of transport through a single impurity Anderson model. The  on-site impurity Hamiltonian is 
# $$
# H_{loc} = U n_\uparrow n_\downarrow + \left(\mu-\frac{U}{2} + \epsilon_d\right)(n_\uparrow + n_\downarrow),
# $$
# where $U$ is the local interaction, and $\epsilon_d$ the bare level energy. The impurity is coupled to two infinite metallic leads, the left ($L$) and right ($R$) bath, at a voltage bias $V$. The corresponding embedding self-energy $\Sigma_{\rm bath}=\Sigma_{L}+\Sigma_{R}$ is defined via a smooth box density of states
# $$
# A_\Sigma(\omega) = \frac{\Gamma/\pi}{(1+e^{\nu(\omega-\omega_c)})(1+e^{-\nu(\omega_c+\omega)})}.
# $$
# $\Gamma$ is the hybridization strength (used as energy unit in the following), $\omega_c$ is the half bandwidth, and $\nu$ a smoothening parameter; $\Sigma^<_{L,R}$ then has distribution functions with different chemical potentials $\mu_{L/R}=\pm V/2$. The noninteracting impurity Green's function $G_{0}$ is therefore determined through the steady-state variant of the Dyson equation $G_0 =( i\partial_t + \mu - \Sigma_{\rm bath})^{-1}$, while the interacting Green's function $G$ is obtained from a Dyson equation  $G =( i\partial_t + \mu - \Sigma_{\rm bath}+\Sigma_U)^{-1}$, with an additional self-energy $\Sigma_U$ due to interactions. In the example below, we approximate  $\Sigma_U$ by a (non self-consistent) 2nd order perturbation theory
# $$
# \Sigma(t,t') = U^2 G_0(t,t') G_0(t,t') G_0(t',t),
# $$
# evaluated in the steady state. 

# ## Files
# 
# - `programs/ness2_siam.cpp`: Source code for the real-time calculation. 
# - `utils/demo_ness2_siam.ipynb`: This notebook.
# - `utils/demo_ness2_siam.py`: Python script, does the same as this Notebook.
# 
# 
# ## Running/Compilation:
# 
# For the installation, follow the instructions for the installation of the `NESSi` examples. The programs must be compiled with `hdf5=ON` and `ness=ON`. Copy the notebook to a working directory, and adapt the path of the exectutables (see below). This script depends on the Python utitilies for reading hdf5 Green's function files. Make sure that the `libcntr/python3` is part of the `PYTHONPATH`, e.g., using the `sys.path.append` code below

# In[ ]:


import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import h5py
import sys


# In[ ]:


# === Global paths ===
# adapt this as nessesary:
path_to_libcntr = "/path/to/libcntr/"  
path_to_examples = "/path/to/nessi/examples/"  
# directory to store input and output [relative to directory from which you run this notebook]
output_folder = "./"
os.makedirs(output_folder, exist_ok=True)

# exectutables --- 
executable = path_to_examples + "exe/ness2_siam.x"

#Import python modules form nessi/
sys.path.append(path_to_libcntr + "python3/") #update PYTHONPATH
from ReadNESS import *


# The following function is used to run the code for given parameters:
# 
# * `U`: Interaction 
# * `beta`: Inverse temperature
# * `Voltage`: Voltage
# * `epsd`: bare level energy
# * `Nft`: Size of the Fourier domain
# * `h`: Time-step
# * `eta`: Possible regularization of the Dyson equation (not needed here).
# 
# The output of this program will be written to an hdf5 file **output_folder/[run_name]_results.h5**

# In[ ]:


def run_program(run_name, h, Nft, beta, epsd, U, Voltage, eta):
    #Writes input file and runs the program, storing results with run_name prefix.    
    # Paths and names
    input_filename = os.path.join(output_folder, run_name + "_input.txt")
    out_file = output_folder + run_name + "_results.h5"    
    log_file = os.path.join(output_folder, run_name + ".log")
    # Write input file
    with open(input_filename, 'w') as f:
        f.write(f"__U= {U}\n")
        f.write(f"__beta= {beta}\n")
        f.write(f"__V= {Voltage}\n")
        f.write(f"__epsd= {epsd}\n")
        f.write(f"__Nft= {Nft}\n")        
        f.write(f"__h= {h}\n")
        f.write(f"__eta= {eta}\n")

    # Run the executable and redirect stdout/stderr to logfile
    with open(log_file, "w") as log:
        try:
            subprocess.run(
                [executable, input_filename, out_file],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,  # merge stderr into stdout
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error during execution, see {log_file} for details.")
            raise


# In[ ]:


# a first run at some given parameters:
h=0.01
Nft=2**14
beta=20
epsd=0
U=6
Voltage=2
eta=0
run_name="run001"
run_program(run_name, h, Nft, beta, epsd, U, Voltage, eta)
out_file = output_folder + run_name + "_results.h5"


# ### Analysis of real-time Green's function and Self-energy
# 
# In order to extract the Green's functions on the FFT time and frequency grid, we use the `read_herm_matrix_ness_h5file` function from `readNESS`:
# 
# `G=read_herm_matrix_ness_h5file(file,key="label")` looks for a Green's function under a group `label/` in the HDF5 file file and returns a class `G` such that `G.ret_t`, `G.les_t`, `G.ret_w`, `G.les_t` are arrays which contain $G^{R|<}(t)$ and $G^{R|<}(\omega)$. The data layout follows the FFT convention: 
# 
# * `G.ret_t[l]` = $G^R(lh)$ for `l=0,...,Nft/2-1` and $G^R((N_{ft}-l)h)$ for `l=Nft/2,...,Nft-1`
# * `G.ret_w[l]` = $G^R(l\Delta_\omega)$ for `l=0,...,Nft/2-1` and $G^R((N_{ft}-l)\Delta_\omega)$ for `l=Nft/2,...,Nft-1`
# 
# (similar for the lesser component) . For plotting, one must rearrange the arrays such that times and frequencies are in increasing order.

# In[ ]:


SL=read_herm_matrix_ness_h5file(out_file, key="SL")
SR=read_herm_matrix_ness_h5file(out_file, key="SR")
G=read_herm_matrix_ness_h5file(out_file, key="G")
G0=read_herm_matrix_ness_h5file(out_file, key="G0")

#freq and time arrays for plotting
dw=2*np.pi/(h*Nft)
freqs = np.arange(-Nft//2, Nft//2)*dw
times = np.arange(-Nft//2, Nft//2)*h
# arap arrays:
sLret=np.concatenate((SL.ret_w[Nft//2:], SL.ret_w[:Nft//2]))
sLles=np.concatenate((SL.les_w[Nft//2:], SL.les_w[:Nft//2]))
sRles=np.concatenate((SR.les_w[Nft//2:], SR.les_w[:Nft//2]))
Gret=np.concatenate((G.ret_w[Nft//2:], G.ret_w[:Nft//2]))
Gles=np.concatenate((G.les_w[Nft//2:], G.les_w[:Nft//2]))
G0ret=np.concatenate((G0.ret_w[Nft//2:], G0.ret_w[:Nft//2]))
G0les=np.concatenate((G0.les_w[Nft//2:], G0.les_w[:Nft//2]))


#plot spectrum and lesser for Sigma, G, G0

#fig, axs = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
fig, axs = plt.subplots(1, 3, figsize=(9, 3.5), sharex=False)
axs.flatten()
nmax=100
# === Top left: Re/Im[G^R(t)]
axs[0].plot(freqs,-2*np.imag(sLret), '-', label=r"$\pi A_\Sigma(\omega)$")
axs[0].plot(freqs,np.imag(sLles), '--', label=r"Im $\Sigma_L^<(\omega)$")
axs[0].plot(freqs,np.imag(sRles), '--', label=r"Im $\Sigma_R^<(\omega)$")

#axs[0].plot(np.imag(SL.les_w/(2*np.pi)), '-', label=r"Im \Sigma_L^<(\omega)")
#axs[0].plot(np.imag(SR.les_w/(2*np.pi)), '-', label=r"Im \Sigma_L^<(\omega)")
axs[0].set_ylabel("")
axs[0].set_xlabel(r"$\omega$")
#axs[0].set_title("")
axs[0].legend(loc="upper right", bbox_to_anchor=(0.93, 0.4))
#axs[0].grid(True, linestyle='-', linewidth=0.5)
#axs[0].set_xlim(0,tmax)
axs[0].set_xlim(-12,12)

# === 
axs[1].plot(freqs,-2*np.imag(Gret), 'b-', label=r"$\pi A(\omega)$")
axs[1].plot(freqs,-2*np.imag(G0ret), 'r-', label=r"$\pi A_0(\omega)$")
axs[1].plot(freqs,np.imag(Gles), 'b--', label=r"$G^<(\omega)$")
axs[1].plot(freqs,np.imag(G0les), 'r--', label=r"$G_0^<(\omega)$")


axs[1].set_ylim(-0.1,1.1)
axs[1].set_xlim(-12,12)
axs[1].set_ylabel("")
axs[1].set_xlabel(r"$\omega$")
axs[1].legend()

# === 
axs[2].plot(freqs,-0.5*np.imag(Gles)/(np.imag(Gret)+1e-8), 'b-', label=r"$F(\omega)$")
axs[2].plot(freqs,-0.5*np.imag(G0les)/(np.imag(G0ret)+1e-8), 'r-', label=r"$F_0(\omega)$")
axs[2].set_ylim(-0.1,1.18)
axs[2].set_xlim(-3,3)
axs[2].set_ylabel("")
axs[2].set_xlabel(r"$\omega$")
axs[2].legend()
#axs[0].grid(True, linestyle='-', linewidth=0.5)
#axs[0].set_xlim(0,tmax)

axs[0].text(0.02, .95, "a)", transform=axs[0].transAxes, fontsize=12, fontweight='bold', va='top')
axs[1].text(0.02, .95, "b)", transform=axs[1].transAxes, fontsize=12, fontweight='bold', va='top')
axs[2].text(0.02, .95, "c)", transform=axs[2].transAxes, fontsize=12, fontweight='bold', va='top')

plt.tight_layout()

#saving figure for paper:
png_file = os.path.join(output_folder, "ness2_siam_spectra.png")
plt.tight_layout()
plt.savefig(png_file, dpi=300, bbox_inches="tight")  # high resolution
plt.show()


# ### Current Voltage characteristics:
# 
# The current, which is defined by the rate of particle transfer form the left to the right reservoir, $J= \frac{dN_L}{dt}= -\frac{dN_R}{dt}$, is given by equal-time convolution 
# $$
# J = 2\left(\Sigma_L\ast G -  G\ast\Sigma_L\right)^<(t,t).
# $$
# This relation is in general not satisfied away from half-filling in bare (non-conserving) second order perturbation theory, which we make use of for calculating the self-energy.
# In the steady state implementation, this can be evaluated using the ``ness2::convolution_density_matrix`` method.

# In[ ]:


# Now we do a scan to get the current-voltage characteristics
Vs = np.linspace(0, 25, 21)

h = 0.01
Nft = 2**14
beta = 20
epsd = 0
eta = 0

U_list = [0, 2, 4, 6]
results = {}
# ====== run the simulations 
for U in U_list:
    Js = np.zeros_like(Vs, dtype=float)
    run_name = "run001"
    for idx, V in enumerate(Vs):
        run_program(run_name, h, Nft, beta, epsd, U, V, eta)
        out_file = os.path.join(output_folder, run_name + "_results.h5")
        with h5py.File(out_file, "r") as f:
            Js[idx] = f["Current"][()][0]  # extract scalar
    results[U] = Js

fig, axs = plt.subplots(1, 1, figsize=(5, 5), sharex=False)

# ====== plot:

for U in U_list:
    axs.plot(Vs, results[U], "o-", lw=1.5, ms=4, label=f"U = {U}")
axs.set_xlabel(r"$V$")
axs.set_ylabel(r"$J$")
#plt.title("I–V characteristics vs interaction U")
#plt.grid(True, ls="--", lw=0.5)
axs.legend(loc="upper right", bbox_to_anchor=(0.93, 0.5))


plt.tight_layout()

#saving figure for paper:
png_file = os.path.join(output_folder, "ness2_siam_curent.png")
plt.savefig(png_file, dpi=300, bbox_inches="tight")  # high resolution

plt.show()


# In[ ]:


## prepare the figure for the manuscript

fig, axs = plt.subplots(2, 2, figsize=(9, 6), sharex=False)
#axs.flatten()
nmax=100
# === Top left: Re/Im[G^R(t)]
axs[0,0].plot(freqs,-2*np.imag(sLret), '-', label=r"$\pi A_\Sigma(\omega)$")
axs[0,0].plot(freqs,np.imag(sLles), '--', label=r"Im $\Sigma_L^<(\omega)$")
axs[0,0].plot(freqs,np.imag(sRles), '--', label=r"Im $\Sigma_R^<(\omega)$")
axs[0,0].set_ylabel("")
axs[0,0].set_xlabel(r"$\omega$")
#axs[0,0].set_title("")
axs[0,0].legend(loc="upper right", bbox_to_anchor=(0.9, 0.45))
#axs[0,0].grid(True, linestyle='-', linewidth=0.5)
#axs[0,0].set_xlim(0,tmax)
axs[0,0].set_xlim(-12,12)

# === 
axs[1,0].plot(freqs,-2*np.imag(Gret), 'b-', label=r"$\pi A(\omega)$")
axs[1,0].plot(freqs,-2*np.imag(G0ret), 'r-', label=r"$\pi A_0(\omega)$")
axs[1,0].plot(freqs,np.imag(Gles), 'b--', label=r"Im$G^<(\omega)$")
axs[1,0].plot(freqs,np.imag(G0les), 'r--', label=r"Im$G_0^<(\omega)$")


axs[1,0].set_ylim(-0.1,1.1)
axs[1,0].set_xlim(-12,12)
axs[1,0].set_ylabel("")
axs[1,0].set_xlabel(r"$\omega$")
axs[1,0].legend()

# === 
axs[1,1].plot(freqs,-0.5*np.imag(Gles)/(np.imag(Gret)+1e-8), 'b-', label=r"$F(\omega)$")
axs[1,1].plot(freqs,-0.5*np.imag(G0les)/(np.imag(G0ret)+1e-8), 'r-', label=r"$F_0(\omega)$")
axs[1,1].set_ylim(-0.1,1.18)
axs[1,1].set_xlim(-3,3)
axs[1,1].set_ylabel("")
axs[1,1].set_xlabel(r"$\omega$")
axs[1,1].legend()
#axs[0].grid(True, linestyle='-', linewidth=0.5)
#axs[0].set_xlim(0,tmax)

axs[0,0].text(0.02, .95, "a)", transform=axs[0,0].transAxes, fontsize=12, fontweight='bold', va='top')
axs[1,0].text(0.02, .95, "c)", transform=axs[1,0].transAxes, fontsize=12, fontweight='bold', va='top')
axs[1,1].text(0.02, .95, "d)", transform=axs[1,1].transAxes, fontsize=12, fontweight='bold', va='top')
axs[0,1].text(0.02, .95, "b)", transform=axs[0,1].transAxes, fontsize=12, fontweight='bold', va='top')

for U in U_list:
    axs[0,1].plot(Vs, results[U], "o-", lw=1.5, ms=4, label=f"U = {U}")
axs[0,1].set_xlabel(r"$V$")
axs[0,1].set_ylabel(r"$J$")
#plt.title("I–V characteristics vs interaction U")
#plt.grid(True, ls="--", lw=0.5)
axs[0,1].legend(loc="upper right", bbox_to_anchor=(0.93, 0.5))

plt.tight_layout()

#saving figure for paper:
png_file = os.path.join(output_folder, "ness2_siam_01.png")
plt.tight_layout()
plt.savefig(png_file, dpi=300, bbox_inches="tight")  # high resolution
plt.show()



# ### Convergence analysis
# 
# * We first fix a small timestep `h=0.001`, and perform simulations with increasing `Nft`. We analyse the error $\epsilon^{R,<}(N_{ft}) = |J_{Nft} - J_{ref}|_2 $, where $J_{Nft}$ is the result of a simulation with given $N_{ft}$, and  the reference $J_{ref}$ is  the result with the largest `Nft=2**20`.
#   
# * We then perform simulations at fixed `Nft` * `h` with different values of `Nft`. These are essentially simulations with varying `h` at fixed largest time  ($t_c$ ~ `Nft/2` * `h`). We choose `Nft` * `h_ness` = `2**20` * `0.001`.

# In[ ]:


# Convergence analysis:
# --- Choose the physical point ---
U, V = 6.0, 2.0
beta = 20.0
epsd = 0.0
eta = 0.0

# --- Convergence sweep settings ---
h_fixed = 0.001                # for sweep (1)
kmax=20
Nft_list_1 = [2**k for k in range(8, kmax+1)]  

T_ref = h_fixed * Nft_list_1[-1]          # fix T = h*Nft for sweep (2)
h_list_2 = [T_ref / (nft) for nft in Nft_list_1]  # decreasing h, Nft grows

def read_current(h5path):
    with h5py.File(h5path, "r") as f:
        return float(np.asarray(f["Current"])[0])

def run_and_get_current(run_name, h, Nft):
    out_file = os.path.join(output_folder, run_name + "_results.h5")
    run_program(run_name, h, Nft, beta, epsd, U, V, eta)
    return read_current(out_file)

# ---- (1) Fixed h, increase Nft ----
J_vals_1 = []
for Nft in Nft_list_1:
    rn = f"conv1_h{h_fixed:g}_N{Nft}"
    J_vals_1.append(run_and_get_current(rn, h_fixed, Nft))
J_vals_1 = np.array(J_vals_1)
J_ref_1 = J_vals_1[-1]                         # take largest Nft as reference
err_1 = np.abs(J_vals_1 - J_ref_1)

# ---- (2) Fixed T=h*Nft, decrease h (increase Nft) ----
h_list_2 = np.array(h_list_2)
Nft_list_2 = Nft_list_1    # ensure integer powers of two from construction
J_vals_2 = []
for h2, Nft2 in zip(h_list_2, Nft_list_2):
    rn = f"conv2_T{T_ref:g}_h{h2:.6g}_N{Nft2}"
    J_vals_2.append(run_and_get_current(rn, h2, Nft2))
J_vals_2 = np.array(J_vals_2)
J_ref_2 = J_vals_2[-1]                         # take smallest h (largest Nft) as reference
err_2 = np.abs(J_vals_2 - J_ref_2)


# In[ ]:


# ---- Plot: two panels, log error ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5))

# Panel 1: fixed h, increase Nft
ax1.plot(Nft_list_1[:-1], err_1[:-1]+1e-15, "o-", lw=1.5)
ax1.set_xscale("log",base=2)
ax1.set_yscale("log")
ax1.set_xlabel(r"$N_{\mathrm{ft}}$")
ax1.set_ylabel(r"$|J - J_{\mathrm{ref}}|$")
#ax1.set_title("Convergence vs $N_{\\rm ft}$ (fixed $h$)")
ax1.grid(True, ls="--", lw=0.5)

# Panel 2: fixed T, decrease h
ax2.plot(h_list_2[:-1], err_2[:-1]+1e-15, "o-", lw=1.5)
ax2.plot(h_list_2[3:-1], h_list_2[3:-1]**2*0.01, "--", lw=1.5)

ax2.set_xscale("log", base=10)
ax2.set_yscale("log")
ax2.invert_xaxis()  # smaller h to the right
ax2.set_xlabel(r"$h$")
ax2.set_ylabel(r"$|J - J_{\mathrm{ref}}|$")
#ax2.set_title("Convergence vs $h$ (fixed $T$)")
ax2.grid(True, ls="--", lw=0.5)

ax1.text(0.1, .8, "a)", transform=ax1.transAxes, fontsize=12, fontweight='bold', va='top')
ax1.text(0.1, .2, r"fixed $h=0.001$", transform=ax1.transAxes, fontsize=12, va='top')
ax2.text(0.1, .8, "b)", transform=ax2.transAxes, fontsize=12, fontweight='bold', va='top')
ax2.text(0.1, .2, r"fixed $N_{ft}h=0.001\times 2^{20}$", transform=ax2.transAxes, fontsize=12, va='top')


#saving figure for paper:
png_file = os.path.join(output_folder, "ness2_siam_convergence.png")
plt.tight_layout()
plt.savefig(png_file, dpi=300, bbox_inches="tight")  # high resolution
plt.show()



# In[ ]:


# clear all data files produced by this worksheet:
os.system(f"rm {output_folder}run001* {output_folder}conv1* {output_folder}conv2*")


# In[ ]:




