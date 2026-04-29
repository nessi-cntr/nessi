#!/usr/bin/env python
# coding: utf-8

# # Demo: DMFT in the steady state

# ## Overview:
# 
# This test program provides an example for the implementation of the steady-state Kadanoff Baym Equations and a benchmark to the real time propagation. We consider a simple Hubbard model on a Bethe lattice 
# $$
# H = -J \sum_{\langle i,j \rangle, \sigma} c_{i \sigma}^\dagger c_{j \sigma} + U \sum_i \Big(n_{i \uparrow}-\tfrac12\Big)\Big(n_{i \downarrow}-\tfrac12\Big) + \mu \sum_i \big(n_{i \uparrow}+n_{i \downarrow}\big).
# $$
# Here, $ c_{i \sigma}^\dagger $ and $ c_{i \sigma} $ are the fermionic creation and annihilation operators at sites $ i,j $ with spin $ \sigma $ and $ n_{i \sigma} = c_{i \sigma}^\dagger c_{i \sigma} $ the number operator. We solve the problem using DMFT, using second-order perturbation theory as an impurity solver. For the Bethe lattice with bandwidth $W=4$, this requires a solution of the Dyson equation equation
# \begin{align}
# G^{-1} &= (i\partial_t + \mu - \Sigma - \Delta),
# \end{align}
# for the local Green's function $G$, with the self-consistency condition $\Delta=G$. The self-energy will be approximated using 2nd order parturbation theory
# $$
# \Sigma(t,t') = U^2 \mathcal{G}(t,t') \mathcal{G}(t,t') \mathcal{G}(t',t),
# $$
# where either $\mathcal{G}=G$ (self-consistent perturbation theory) or $\mathcal{G}=(i\partial_t + \mu -\Delta)^{-1}$ (iterated perturbation theory, IPT).
# 
# We solve the problem in thermal equilibrium (temperature $T=1/\beta$), in order to compare the implemenation using the real-time propagation and the direct steady-state approach.
# 

# ## Files
# 
# - `programs/ness2_bethe_prop.cpp`: Source code for the real-time calculation. 
# - `programs/ness2_bethe.cpp`: Source code for the steady-state calculation. 
# - `utils/demo_ness2_bethe.ipynb`: This notebook.
# - `utils/demo_ness2_bethe.py`: Python script, does the same as this Notebook.
# 
# 
# ## Running/Compilation:
# 
# * For the installation, follow the instructions for the installation of the `NESSi` examples. The programs must be compiled with `hdf5=ON` and `ness=ON`
# * Copy the notebook to a working directory, and adapt the path of the exectutables (see below).
# * This script depends on the Python utitilies for reading hdf5 Green's function files. Make sure that the `libcntr/python3` is part of the `PYTHONPATH`, e.g., using the `sys.path.append` code below

# In[ ]:


# All python modules
import numpy as np
import os
import subprocess
import h5py
import sys
# === for plotting
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec


# In[ ]:


# === Path ===  adapt this section as nessesary:
path_to_libcntr = "/path/to/nessi/libcntr/"  
path_to_examples = "/path/to/nessi/examples/"
# directory to store ionput and output [relative to directory from which you run this notebook]
output_folder = "./"
os.makedirs(output_folder, exist_ok=True)

# exectutables --- 
executable_ness = path_to_examples + "exe/ness2_bethe.x"
executable_cntr = path_to_examples + "exe/ness2_bethe_prop.x"
executable_trunc = path_to_examples + "exe/ness2_bethe_trunc.x"

#Import python modules form nessi/
sys.path.append(path_to_libcntr + "python3/") #update PYTHONPATH
from ReadNESS import *
from ReadCNTR import *
from ReadCNTRhdf5 import *


# ## Real time simulation: 
# 
# The self-consistent KBEs are solved using the standard `NESSi` `dyson` routines, with a accuracy `SolveOrder=MAX_SOLVE_ORDER=5`
# 
# **Parameters:**
# 
# * `U`: Interaction 
# * `mu`: Chemical potential
# * `beta`: Inverse temperature
# * `nt_cntr`: Number of timesteps in the full KBE simulation.
# * `h_cntr`: Timestep
# * `ntau`: Number of timesteps on the imaginary axis
# * `BootstrapMaxIter`: Max number of iteratins for the bootstrapping phase (time 1,...,SolveOrder)
# * `BootstrapMaxErr`: Error cutoff for the bootstrapping phase
# * `CorrectorSteps`: Fixed Number of selfconsistency iterations at each timstep
# * `ipt_flag`: Whether to use IPT selfenergy (`ipt_flag=1`) or not.
# 
# 
# (The parameters $U$, $\mu$, $\beta$, `ipt_flag` will be the same for the NESS simulation below.)
# 
# The output of this program will be written to an hdf5 file **output_folder/[run_name]_results_cntr.h5**

# In[ ]:


def run_program_cntr(run_name,nt, h, beta, U, mu, ipt_flag, ntau, BootstrapMaxIter, BootstrapMaxErr, CorrectorSteps):  
    input_file = os.path.join(output_folder, f"{run_name}_input_cntr.txt")
    out_file   = os.path.join(output_folder, f"{run_name}_results_cntr.h5")
    log_file   = os.path.join(output_folder, f"{run_name}_cntr.log")

    # write input file (same key names/order as your example)
    with open(input_file, 'w') as f:
        f.write(f"__nt={nt}\n")
        f.write(f"__h={h}\n")
        f.write(f"__beta={beta}\n")
        f.write(f"__U={U}\n")
        f.write(f"__mu={mu}\n")
        f.write(f"__ipt_flag={ipt_flag}\n")
        f.write(f"__ntau={ntau}\n")
        f.write(f"__BootstrapMaxIter={BootstrapMaxIter}\n")
        f.write(f"__BootstrapMaxErr={BootstrapMaxErr}\n")
        f.write(f"__CorrectorSteps={CorrectorSteps}\n")

    # run executable, redirecting all output to a logfile
    with open(log_file, "w") as log:
        try:
            subprocess.run(
                [executable_cntr, input_file, out_file],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True
            )
        except subprocess.CalledProcessError:
            print(f"CNTR run failed. See log: {log_file}")
            raise


# In[ ]:


# === the following params will be the same for NESS simulation
run_name="run001"
U=4
mu=0 # half filling
beta=10
ipt_flag=1
h_cntr=0.02  # rather small for the real-time simulation, to make sure results are accurate
nt_cntr=1000 # real-time solved on time-steps [0,...,nt_cntr]
ntau=1000    # dtau=beta/ntau ~ h is a reasonable choice
BootstrapMaxIter=100 
BootstrapMaxErr=1e-8
CorrectorSteps=4
# ==========================
# running the real-time program; for (nt=1000,ntau=1000) this may take about a minute
run_program_cntr(run_name, nt_cntr, h_cntr, beta, U, mu, ipt_flag, ntau, BootstrapMaxIter, BootstrapMaxErr, CorrectorSteps)
out_file_cntr = output_folder + run_name + "_results_cntr.h5"


# In[ ]:


# Quick check: Extract DMFT convergence error, to confirm convergence:

errs = []
with h5py.File(out_file_cntr, "r") as f:
    itermax = int(f["itermax"][0])
    #print(f"Converged after = {itermax} iterations")
    for it in range(itermax + 1):
        errs.append(float(f[f"it{it}/err"][0])) 

# Plot
fig, axs = plt.subplots(1, 1, figsize=(5, 4), sharex=False)

its = np.arange(itermax + 1)
axs.plot(its, errs, 'k-', linewidth=1.0)
axs.scatter(its, errs, s=30)
axs.set_yscale('log')
axs.set_xlabel("iteration")
axs.set_ylabel("error")
axs.set_title("Convergence")

plt.show()


# ### Analysis of real-time Green's function and Self-energy
# 
# The following code extracts the Green's functions and self-energy such that 
# 
# * `Gret_cntr[l]` = $G^R(t_{max}, t_{max} -t_l)$, `Gles_cntr[l]` = $G^<(t_{max}, t_{max} -t_l)$,  
# * `Sret_cntr[l]` = $\Sigma^R(t_{max}, t_{max} -t_l)$,  `Sret_cntr[l]` = $\Sigma^R(t_{max}, t_{max} -t_l)$,
# 
# at the last timeslice $t_{max}=$ `nt_cntr` * `h_cntr` with $t_l= l$*`h_cntr`

# In[ ]:


def read_last_slice(file, key):
    with h5py.File(file, 'r') as fd:
        G = read_gf_group(fd[key])   # from ReadCNTRhdf5
        ret = G.ret[:, :, 0, 0]      # scalar case
        les = G.les[:, :, 0, 0]

        # Fill lower triangle of lesser via Hermitian relation (t1>t2)
        i_low = np.tril_indices_from(les, k=-1)
        les[i_low] = -les.T[i_low].conj()

        nt = ret.shape[0] - 1                # last physical time index
        ret_slice = np.flip(ret[nt, :nt+1])   # t2 = 0..nt-1, reversed
        les_slice = np.flip(les[nt, :nt+1])
        return ret_slice, les_slice

Sret_cntr, Sles_cntr = read_last_slice(out_file_cntr, "S")
Gret_cntr, Gles_cntr = read_last_slice(out_file_cntr ,"G")


# In[ ]:


# Plot GF and sigma on the last timeslice
times_cntr=np.arange(nt_cntr+1) * h_cntr

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes = axes.flatten()

axes[0].plot(times_cntr,np.imag(Gret_cntr), 'b-',label=r'Im$G^R(t_{max},t_{max}-s)$')
axes[0].plot(times_cntr,np.real(Gret_cntr), 'b--', label=r'Re$G^R(t_{max},t_{max}-s)$')
axes[0].plot(times_cntr,np.imag(Gles_cntr), 'r-',label=r'Im$G^<(t_{max},t_{max}-s)$')
axes[0].plot(times_cntr,np.real(Gles_cntr), 'r--', label=r'Re$G^<(t_{max},t_{max}-s)$')
axes[0].set_xlabel(r"$s$")
axes[0].set_ylabel(r"")
axes[0].legend()

axes[1].plot(times_cntr,np.imag(Sret_cntr), 'b-',label=r'Im$\Sigma^R(t_{max},t_{max}-s)$')
axes[1].plot(times_cntr,np.real(Sret_cntr), 'b--', label=r'Re$\Sigma^R(t_{max},t_{max}-s)$')
axes[1].plot(times_cntr,np.imag(Sles_cntr), 'r-',label=r'Im$\Sigma^<(t_{max},t_{max}-s)$')
axes[1].plot(times_cntr,np.real(Sles_cntr), 'r--', label=r'Re$\Sigma^<(t_{max},t_{max}-s)$')
axes[1].set_xlabel(r"$s$")
axes[1].set_ylabel(r"")
axes[1].legend()

#plt.tight_layout()
plt.show()


# ## NESS simulation
# 
# The steady-state KBEs are solved using the standard steady-state Dyson routine `ness2::dyson`.
# 
# **Parameters:**
# 
# (The parameters $U$, $\mu$, $\beta$, `ipt_flag` are the same as for the cntr simulation above.
# 
# * `Nft`: Size of the Fourier domain
# * `h_ness`: timestep
# * `N_it`: Number of self-consistent DMFT iterations
# * `mixing`: Linear mixing for DMFT convergence (more relevant than for the real-time simulation)
# * `out_every`: Flag to print full Green's function at every given DMFT iteration
# 
# **Note:** The maximal time represented by the grid is $t_{cut}$ = `h_ness` * (`Nft/2-1`). This should be sufficiently large such that steady-state Green's functions $G^{R,<}(t)$ are close to zero for relative times $|t|>t_{cut}$. A systematic convergence check will be done below.
# 

# In[ ]:


def run_program_ness(run_name, U, beta, Nft, mu,ipt_flag, h, Niter, errmax,  mixing, out_every):
    #Writes input file and runs the program, storing results with run_name prefix.    
    # Paths and names
    input_filename = os.path.join(output_folder, run_name + "_input.txt")
    out_file = output_folder + run_name + "_results_ness.h5"    
    log_file = os.path.join(output_folder, run_name + ".log")
    # Write input file
    with open(input_filename, 'w') as f:
        f.write(f"__U= {U}\n")
        f.write(f"__beta= {beta}\n")
        f.write(f"__Nft= {Nft}\n")        
        f.write(f"__mu= {mu}\n")
        f.write(f"__ipt_flag= {ipt_flag}\n")
        f.write(f"__h= {h}\n")
        f.write(f"__Niter= {Niter}\n")
        f.write(f"__errmax= {errmax}\n")
        f.write(f"__mixing= {mixing}\n")
        f.write(f"__out_every= {out_every}\n")
    # Run the executable and redirect stdout/stderr to logfile
    with open(log_file, "w") as log:
        try:
            subprocess.run(
                [executable_ness, input_filename, out_file],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,  # merge stderr into stdout
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error during execution, see {log_file} for details.")
            raise


# In[ ]:


# A first run at some parameters (systematoic convergece checks will be cone below)
Nft=2**13     # convergence will be checked below
h_ness=h_cntr # convergence will be checked below
Niter=100
mixing=0.6
errmax=1e-7
out_every=1

run_program_ness(run_name, U, beta, Nft, mu,ipt_flag, h_ness, Niter, errmax,  mixing, out_every)
out_file_ness = output_folder + run_name + "_results_ness.h5"


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

# ### Analysis of the DMFT convergence
# 
# Plot the DMFT spectrum for differenet iterations. When playing with the `mixing` parameter, one would see that a small mixing helps convergence significantly, especially for late times, where the Green's function is very small.
#   

# In[ ]:


# read data such that Gs[it] is Green's function at DMFT iteration it:
errs, Gs = [],[]
with h5py.File(out_file_ness, "r") as f:
    itermax = int(f["itermax"][0])
    #print(f"Converged after = {itermax} iterations")
    for it in range(itermax + 1):
        Gs.append(read_herm_matrix_ness_h5file(out_file_ness, key=f"it{it}/G"))    
        errs.append(float(f[f"it{it}/err"][0]))  # scalar

# ==== figure =====

cmap = matplotlib.colormaps['bwr']  # blue → white → red
norm = mcolors.Normalize(vmin=0, vmax=itermax)

#fig, axs = plt.subplots(1, 2, figsize=(9, 4), sharex=False, constrained_layout=True)
fig, axs = plt.subplots(1, 2, figsize=(9, 4), sharex=False)
fig.subplots_adjust(wspace=0.22)

freqs=np.arange(-Nft//2,Nft//2)*(2*np.pi/(Nft*h_ness))

# --- Panel 1: Im G^R ---
for it in range(itermax + 1):
    Gret_ness = np.concatenate((Gs[it].ret_w[Nft // 2:], Gs[it].ret_w[:Nft // 2]))
    axs[0].plot(freqs, -2 * np.imag(Gret_ness), '-', color=cmap(norm(it)))
axs[0].set_ylabel(r"Im$\,G^R(\omega)$")
axs[0].set_xlabel(r"$\omega$")
axs[0].set_xlim(-5, 5)

# --- Panel 2: Im G^< ---
for it in range(itermax + 1):
    Gles_ness = np.concatenate((Gs[it].les_w[Nft // 2:], Gs[it].les_w[:Nft // 2]))
    axs[1].plot(freqs, np.imag(Gles_ness), '-', color=cmap(norm(it)))
axs[1].set_ylabel(r"Im$\,G^<(\omega)$")
axs[1].set_xlabel(r"$\omega$")
axs[1].set_xlim(-5, 1)

# --- Inset: error vs iteration ---
ax_inset = axs[1].inset_axes([0.2,0.65,0.35,0.3])

its = np.arange(itermax + 1)
colors = cmap(norm(its[1:]))
ax_inset.plot(its[1:], errs[1:], '-', linewidth=1.0, color='0.6')
ax_inset.scatter(its[1:], errs[1:], c=colors, s=20)
ax_inset.set_yscale('log')
ax_inset.set_xlabel("iteration", fontsize=8)
ax_inset.set_ylabel("", fontsize=8)
ax_inset.tick_params(axis='both', which='major', labelsize=8)

axs[0].text(0.02, .95, "a)", transform=axs[0].transAxes, fontsize=12, fontweight='bold', va='top')
axs[1].text(0.02, .95, "b)", transform=axs[1].transAxes, fontsize=12, fontweight='bold', va='top')
ax_inset.text(0.5,0.9, "err", transform=ax_inset.transAxes, fontsize=9, va='top')

# --- colorbar ---
cax=fig.add_axes([0.4,0.55,0.01,0.3])
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cax, orientation='vertical', fraction=0.02, pad=0.02)
cbar.set_label("Iteration")

# ===== saving figure for paper ===== 
png_file = os.path.join(output_folder, "ness_bethe_convergence.png")
#plt.tight_layout()
plt.savefig(png_file, dpi=300, bbox_inches="tight")

plt.show()


# ### Analysis of the Green's function
# 
# The following code  plots the congerved  $G^{<|R}(t)$ and $\Sigma^{<|R}(t)$, comparing the real-time result (black dashed lines, on top of the NESS results within the linewidth of the plots). On the log-plot, one can see that the Green's function decays well in the given simulation interval, which is needed for a convergence of the FFT-based Dyson equation.

# In[ ]:


#Here we analyze the converged spectra and plot the ness and cntr results as function of time

# read G and Sigma for converged DMFT solution:
Gness=read_herm_matrix_ness_h5file(out_file_ness, key="G")
Sness=read_herm_matrix_ness_h5file(out_file_ness, key="S")

# Using positive times only for plotting:
times_ness = np.arange(Nft // 2) * h_ness
Gret_ness=Gness.ret_t[:Nft//2]
Gles_ness=Gness.les_t[:Nft//2]
Sret_ness=Sness.ret_t[:Nft//2]
Sles_ness=Sness.les_t[:Nft//2]

# Plot GF and sigma on a slice
fig, axes = plt.subplots(2, 2, figsize=(9, 5.5))
axes = axes.flatten()

axes[0].plot(times_ness,np.imag(Gret_ness), 'b-',label=r'Im$G^R(t_1-t_2)$')
axes[0].plot(times_ness,np.real(Gret_ness), 'b--', label=r'Re$G^R(t_1-t_2)$')
axes[0].plot(times_ness,np.imag(Gles_ness), 'r-',label=r'Im$G^<(t_1-t_2)$')
axes[0].plot(times_ness,np.real(Gles_ness), 'r--', label=r'Re$G^<(t_1-t_2)$')

axes[0].plot(times_cntr,np.imag(Gret_cntr), 'k--')
axes[0].plot(times_cntr,np.real(Gret_cntr), 'k--')
axes[0].plot(times_cntr,np.imag(Gles_cntr), 'k--')
axes[0].plot(times_cntr,np.real(Gles_cntr), 'k--')

axes[0].set_xlabel(r"$t_1-t_2$")
axes[0].set_ylabel(r"")
axes[0].legend()
axes[0].set_xlim(0,14)
#axes[0].set_ylim(-0.1,0.1)

axes[1].plot(times_ness,np.imag(Sret_ness), 'b-',label=r'Im$\Sigma^R(t_1-t_2)$')
axes[1].plot(times_ness,np.real(Sret_ness), 'b--', label=r'Re$\Sigma^R(t_1-t_2)$')
axes[1].plot(times_ness,np.imag(Sles_ness), 'r-',label=r'Im$\Sigma^<(t_1-t_2)$')
axes[1].plot(times_ness,np.real(Sles_ness), 'r--', label=r'Re$\Sigma^<(t_1-t_2)$')

axes[1].plot(times_cntr,np.imag(Sret_cntr), 'k--')
axes[1].plot(times_cntr,np.real(Sret_cntr), 'k--')
axes[1].plot(times_cntr,np.imag(Sles_cntr), 'k--')
axes[1].plot(times_cntr,np.real(Sles_cntr), 'k--')


axes[1].set_xlabel(r"$t_1-t_2$")
axes[1].set_ylabel(r"")
axes[1].legend()
axes[1].set_xlim(0,14)

axes[2].plot(times_ness,np.abs(Gret_ness), 'b-',label=r'$|G^R(t_1-t_2)|$')
axes[2].plot(times_ness,np.abs(Gles_ness), 'r-',label=r'$|G^<(t_1-t_2)|$')

axes[2].plot(times_cntr,np.abs(Gret_cntr), 'k--')
axes[2].plot(times_cntr,np.abs(Gles_cntr), 'k--')

axes[2].set_xlabel(r"$t_1-t_2$")
axes[2].set_ylabel(r"")
axes[2].legend()
#axes[2].set_xlim(0,20)
axes[2].set_yscale('log')

axes[3].plot(times_ness,np.abs(Sret_ness), 'b-',label=r'$|\Sigma^R(t_1-t_2)|$')
axes[3].plot(times_ness,np.abs(Sles_ness), 'r-',label=r'$|\Sigma^<(t_1-t_2)|$')

axes[3].plot(times_cntr,np.abs(Sret_cntr), 'k--')
axes[3].plot(times_cntr,np.abs(Sles_cntr), 'k--')


axes[3].set_xlabel(r"$t_1-t_2$")
axes[3].set_ylabel(r"")
axes[3].legend()
axes[3].set_yscale('log')
#axes[3].set_xlim(0,20)

axes[0].text(0.15, .95, "a)", transform=axes[0].transAxes, fontsize=12, fontweight='bold', va='top')
axes[1].text(0.15, .95, "b)", transform=axes[1].transAxes, fontsize=12, fontweight='bold', va='top')
axes[2].text(0.15, .95, "c)", transform=axes[2].transAxes, fontsize=12, fontweight='bold', va='top')
axes[3].text(0.15, .95, "d)", transform=axes[3].transAxes, fontsize=12, fontweight='bold', va='top')

#saving figure for paper:
png_file = os.path.join(output_folder, "ness_bethe_time_comparison.png")
plt.tight_layout()
plt.savefig(png_file, dpi=300, bbox_inches="tight")  # high resolution


plt.tight_layout()
plt.show()


# ### Convergence analysis
# 
# * We first fix the timestep `h_ness`, and perform simulations with increasing `Nft`. We analyse the error $\epsilon^{R,<}(N_{ft}) = |G^{R,<}_{Nft} - G^{R,<}_{ref}|_2 $, where $G_{Nft}$ is the result of a simulation with given $N_{ft}$, and  the reference $G_{ref}$ is  the result with the largest `Nft`. $\epsilon^{R,<}$ is computed as the $L_2$ norm difference between the real-time data $G^{<|R}(t)$ on the common grid of all functions (i.e, the grid with the smallest `Nft`).
# 
# * We then perform simulations at fixed `Nft` * `h_ness` with different values of `Nft`. These are essentially simulations with varying `h_ness` at fixed largest time  ($t_c$ ~ `Nft/2` * `h_ness`). We choose `Nft` * `h_ness` = `2**kref` * `h_cntr`, such that the time-grid of the steady-state simulation coincides with the grid of the real-time simulation for `Nft`=`2**kref`, while smaller (larger) `Nft` correspond to larger (smaller) `h_ness`.  Choosing `Nft` as powers of $2$, all simulations share a common time-grid, on which the result can be compared. Assuming that the real-time simulation provides an accurate benchmark, we can access the convergence of the steady-state simulation with decreasing `h_ness`.

# In[ ]:


# generate data for convergence analysis at fixed h_ness:
h_ness=h_cntr
kmax=20
Nft_list= [2**k for k in range(8, kmax+1)]  

Nft_min=Nft_list[0]
Nft_max=Nft_list[-1]

Gret_list=[]
Gles_list=[]

for Nft in Nft_list:
    run_program_ness(run_name, U, beta, Nft, mu,ipt_flag, h_ness, Niter, errmax,  mixing, 0)
    print(f"... simulation with Nft={Nft}")
    out_file_ness = output_folder + run_name + "_results_ness.h5"
    Gness=read_herm_matrix_ness_h5file(out_file_ness, key="G")
    Gret_list.append(Gness.ret_t[:Nft_min//2])
    Gles_list.append(Gness.les_t[:Nft_min//2])    


# In[ ]:


# generate data for convergence analysis at fixed h_ness*Nft:

kmax=20
# Nft_ref=2**kref, for Nft_ref, we will have h_ness=h_cntr
# tc ~ Nft_ref//2 *h_cntr must be large anough such that G(t) is essentually zero for |t|>tc 
# (this is confirmed by previous convergence analysis at fixed h_ness)
kref=16 
Nft_list_2= [2**k for k in range(11, kmax+1)] 
h_list_2= [ (h_cntr*(2**kref))/(2**k) for k in range(11, kmax+1)] 

Gret_list_2=[]
Gles_list_2=[]

for idx,Nft in enumerate(Nft_list_2):
    h=h_list_2[idx]
    run_program_ness(run_name, U, beta, Nft, mu,ipt_flag, h, Niter, errmax,  mixing, 0)
    print(f"... simulation with Nft={Nft}")
    out_file_ness = output_folder + run_name + "_results_ness.h5"
    Gness=read_herm_matrix_ness_h5file(out_file_ness, key="G")
    #downsampling to common grid with ,argest timestep:
    fac=Nft//Nft_list_2[0]    
    Gret_list_2.append(Gness.ret_t[:Nft//2:fac])
    Gles_list_2.append(Gness.les_t[:Nft//2:fac])    


# In[ ]:


# ===================
# extract errors for convergence test at fixed h_ness
err_les=np.zeros_like(Nft_list[:-1],dtype='float')
err_ret=np.zeros_like(Nft_list[:-1],dtype='float')
for idx,Nft in enumerate(Nft_list[:-1]):
    delta=np.sqrt(np.sum(np.abs(Gret_list[idx]-Gret_list[-1])**2))
    err_ret[idx]=delta.copy()
    delta=np.sqrt(np.sum(np.abs(Gles_list[idx]-Gles_list[-1])**2))
    err_les[idx]=delta.copy()

# ===================
# extract errors for convergence test at fixed h_ness*Nft
#downsampling G_cntr to the grid with largest timestep
faccntr=2**kref//Nft_list_2[0]    
Gles_cntr_2=Gles_cntr[::faccntr]
Gret_cntr_2=Gret_cntr[::faccntr]
ntlen=len(Gret_cntr_2)
err_les_2=np.zeros_like(Nft_list_2,dtype='float')
err_ret_2=np.zeros_like(Nft_list_2,dtype='float')
for idx,Nft in enumerate(Nft_list_2):    
    delta=np.sqrt(np.sum(np.abs(Gret_list_2[idx][:ntlen]-Gret_cntr_2)**2))
    err_ret_2[idx]=delta.copy()
    delta=np.sqrt(np.sum(np.abs(Gles_list_2[idx][:ntlen]-Gles_cntr_2)**2))
    err_les_2[idx]=delta.copy()


# ===============================================

fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
axes = axes.flatten()

axes[0].plot(Nft_list[:-1],err_ret, 'bo-',label=r'$\epsilon(G^R)$')
axes[0].plot(Nft_list[:-1],err_les, 'ro-',label=r'$\epsilon(G^<)$')

axes[0].set_xlabel(r"$N_{ft}$")
axes[0].set_ylabel(r"$\epsilon^{R},\epsilon^{<}$")
axes[0].legend()
axes[0].set_yscale('log')
axes[0].set_xscale('log',base=2)
#axes[0].set_ylim(-0.1,0.1)
axes[0].axvline(2**kref, color='gray', linestyle='--')


axes[1].plot(h_list_2,err_ret_2, 'bo-',label=r'$\epsilon(G^R)$')
axes[1].plot(h_list_2,err_les_2, 'ro-',label=r'$\epsilon(G^<)$')
axes[1].plot(h_list_2[2:],0.2*np.array(h_list_2[2:])**2, color='gray', linestyle='--')

axes[1].set_xlabel(r"$h_{ness}$")
axes[1].set_ylabel(r"$\epsilon^{R},\epsilon^{<}$")

axes[1].legend()
#axes[0].set_xlim(0,14)
#axes[0].set_ylim(-0.1,0.1)
axes[1].set_yscale('log')
axes[1].set_xscale('log',base=10)
axes[1].invert_xaxis()  # smaller h to the right
#axes[1].axvline(h_ness, color='gray', linestyle='--')

axes[0].text(0.25, .95, "a)", transform=axes[0].transAxes, fontsize=12, fontweight='bold', va='top')
axes[1].text(0.25, .95, "b)", transform=axes[1].transAxes, fontsize=12, fontweight='bold', va='top')

#saving figure for paper:
png_file = os.path.join(output_folder, "ness_bethe_benchmark.png")
plt.tight_layout()
plt.savefig(png_file, dpi=300, bbox_inches="tight")  # high resolution


#plt.tight_layout()
plt.show()


# ## Interfacing NESS with Truncated KBE
# 
# The memory-truncated KBE and steady-state NESS Green's functions are both fully defined by their real-time components $G^R$ and $G^<$, making direct data exchange straightforward. In particular, we can use the steady-state simulation to initialize the memory-truncated Green's functions, to perform a real-time simulation starting from an equilibrium state. This setup allows skipping the imaginary-time initialization step, using NESS output directly to seed the truncated KBE evolution. The workflow is:
# 
# 1. **NESS stage**: Run `ness2_bethe.x` to generate a steady-state solution with timestep $h_{\rm ness}$ and frequency-domain size $N_{\rm ft}$.
# 2. **Truncated KBE stage**: Run `ness2_bethe_trunc.x` to read the NESS output and initialize the memory-truncated Green's functions and self-energies. Further propagate the Green's function using the memory trancated KBE.
# 
# Here, we futher propagate with a time-independent Hamiltonian, such that the solution should remain translationally invariant in time. This can be used to access the accuracy.
# 
# 
# **Note:** Because the error in the NESS simulation decreases more slowly  with the timestep $h$ than in the real-time KBE implementation, it can be beneficial to perform the NESS simulation with a smaller timestep than the one used in the truncated evolution. We will set $h_{\rm cntr}=d\cdot h_{\rm ness}$, with a integer factor $d$.
# 
# **Parameters:**
# 
# - `downsampling` - the factor $d$ defined above  
# - `tc` – memory cutoff  
# - `tmax` – propagation length (timesteps)  
# - `out_every` – output frequency  
# - `CorrectorSteps` – DMFT iterations per timestep  
# 
# Further parameters ($U$, $\beta$, $h_{\rm ness}$, $N_{\rm ft}$, `ipt_flag`) are read directly from the NESS output HDF5 file.
# 

# In[ ]:


def run_program_trunc(run_name, downsampling, tc, tmax, out_every, CorrectorSteps):
    input_filename = os.path.join(output_folder, run_name + "_input.txt")
    ness_file = output_folder + run_name + "_results_ness.h5"    
    out_file = output_folder + run_name + "_results_trunc.h5"    
    log_file = os.path.join(output_folder, run_name + ".log")
    # Write input file
    with open(input_filename, 'w') as f:
        f.write(f"__out_every= {out_every}\n")
        f.write(f"__downsampling= {downsampling}\n")
        f.write(f"__tc= {tc}\n")
        f.write(f"__tmax= {tmax}\n")
        f.write(f"__CorrectorSteps= {CorrectorSteps}\n")

    # Run the executable and redirect stdout/stderr to logfile
    with open(log_file, "w") as log:
        try:
            subprocess.run(
                [executable_trunc, input_filename, ness_file, out_file],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,  # merge stderr into stdout
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error during execution, see {log_file} for details.")
            raise

# Extract the error as difference  err^{R,<}(t) = max_s (|G^{R,<}(t,t-s)-G^{R,<}(t_0,t_0-s)|) 
# between G at timeslice t and initial timslice t0=0 (which bhas been set by the NESS initialization)
def evolve_and_get_error(run_name, downsampling, tc, tmax, out_every, CorrectorSteps):
    out_file_trunc = output_folder + run_name + "_results_trunc.h5"
    Gness=read_herm_matrix_ness_h5file(out_file_trunc, key=f"t0/G")
    gret0=Gness.ret_t[:tc+1]
    gles0=Gness.les_t[:tc+1]
    tslices = np.arange(0,tmax+1,out_every)
    err_ret=np.zeros_like(tslices,dtype='float')
    err_les=np.zeros_like(tslices,dtype='float')
    for idx,t in enumerate(tslices):
        Gness=read_herm_matrix_ness_h5file(out_file_trunc, key=f"t{t}/G")
        err_ret[idx]=np.max(np.abs(gret0-Gness.ret_t[:tc+1]))#np.sqrt(np.sum(np.abs(gret0-Gness.ret_t[:tc+1])**2))
        err_les[idx]=np.max(np.abs(gles0-Gness.les_t[:tc+1]))#np.sqrt(np.sum(np.abs(gles0-Gness.les_t[:tc+1])**2))
    return tslices,err_ret,err_les


# In[ ]:


CorrectorSteps=4
tmax=10000
out_every=50
U=1
beta=2
ipt_flag=0

# NOTE: this block will take longer, due to the long real-time simulations


##### ------ 1st simulation: 

run_name="trunc01"
h_cntr=0.04
tc=1000
downsampling=1
Nft=2**15
h_ness=h_cntr/downsampling
run_program_ness(run_name, U, beta, Nft, mu,ipt_flag, h_ness, Niter, errmax,  mixing, 0)
print(f"finished NESS simulation with Nft={Nft} and h={h_ness}")
print(f"starting truncated evolution at tc={tc} using h={h_cntr}")
run_program_trunc(run_name, downsampling, tc, tmax, out_every, CorrectorSteps)
tslices1,err_ret1,err_les1 = evolve_and_get_error(run_name, downsampling, tc, tmax, out_every, CorrectorSteps)

##### ------ 1st simulation: 

run_name="trunc02"
h_cntr=0.04
tc=1000
downsampling=8
Nft=2**18
h_ness=h_cntr/downsampling
run_program_ness(run_name, U, beta, Nft, mu,ipt_flag, h_ness, Niter, errmax,  mixing, 0)
print(f"finished NESS simulation with Nft={Nft} and h={h_ness}")
print(f"starting truncated evolution at tc={tc} using h={h_cntr}")
run_program_trunc(run_name, downsampling, tc, tmax, out_every, CorrectorSteps)
tslices2,err_ret2,err_les2 = evolve_and_get_error(run_name, downsampling, tc, tmax, out_every, CorrectorSteps)


# In[ ]:


fig, axs = plt.subplots(1, 2, figsize=(9, 4), sharex=False)

# =========================
# Plot err^{R,<}(t) = max_s (|G^{R,<}(t,t-s)-G^{R,<}(t_0,t_0-s)|)  as function of time 
h_cntr=0.04
axs[1].plot(tslices1*h_cntr, err_ret1, 'b-', label=r'$h_{\rm ness}=0.04$')
axs[1].plot(tslices1*h_cntr, err_les1, 'b--')

axs[1].plot(tslices2*h_cntr, err_ret2, 'y-', label=r'$h_{\rm ness}=0.005$')
axs[1].plot(tslices2*h_cntr, err_les2, 'y--')


axs[1].set_yscale('log')
axs[1].set_xlabel(r"$t$")
#axs[1].set_ylabel(r'max$_s\,|G(t,t-s)-G(t_0,t_0-s)|$')
axs[1].set_ylabel(r'$\epsilon^R(t)$, $\epsilon^<(t)$')
axs[1].legend()


# =========================
# Plot (|G^{R,<}(t,t-s)-G^{R,<}(t_0,t_0-s)|)  exemplarily for some simesteps t 

#for 3rd
run_name="trunc02"
tc=1000
out_file_trunc = output_folder + run_name + "_results_trunc.h5"
Gness=read_herm_matrix_ness_h5file(out_file_trunc, key=f"t0/G")
gret0=Gness.ret_t[:tc+1]
gles0=Gness.les_t[:tc+1]

times=np.arange(tc+1)*h_cntr
cmap = matplotlib.colormaps['bwr']  # blue → white → red
norm = mcolors.Normalize(vmin=0, vmax=tmax*h_cntr)

for t in [50,10000]: #np.arange(0,tmax+1,10000):
    #print(f"reading t{t}/G")
    Gness=read_herm_matrix_ness_h5file(out_file_trunc, key=f"t{t}/G")
    gret=Gness.ret_t[:tc+1]-gret0
    gles=Gness.les_t[:tc+1]-gles0
    axs[0].plot(times,np.abs(gret), '-',label=f't={t*h_cntr}', color=cmap(norm(t*h_cntr)))
    axs[0].plot(times,np.abs(gles), '--', color=cmap(norm(t*h_cntr)))
axs[0].legend()
axs[0].set_ylabel(r'$|G(t,t-s)-G(t_0,t_0-s)|$')
axs[0].set_xlabel(r'$s$')
axs[0].set_yscale('log')


axs[0].text(0.2, .95, "a)", transform=axs[0].transAxes, fontsize=12, fontweight='bold', va='top')
axs[0].text(0.3, .95, r"$h_{\rm ness}=0.005$", transform=axs[0].transAxes, fontsize=12, va='top')
axs[1].text(0.1, .8, "b)", transform=axs[1].transAxes, fontsize=12, fontweight='bold', va='top')

# --- Shared colorbar ---
#cax=fig.add_axes([0.14,0.15,0.01,0.3])
#sm = cm.ScalarMappable(cmap=cmap, norm=norm)
#sm.set_array([])
#cbar = fig.colorbar(sm, cax=cax, orientation='vertical', fraction=0.02, pad=0.02)
#cbar.set_label(r"timeslice $t$")
plt.tight_layout()

#saving figure for paper:
png_file = os.path.join(output_folder, "ness_bethe_trunc.png")
plt.savefig(png_file, dpi=300, bbox_inches="tight")

plt.show()


# In[ ]:


# clear all data files produced by this worksheet:
os.system(f"rm {output_folder}run001* {output_folder}trunc01* {output_folder}trunc02*")

