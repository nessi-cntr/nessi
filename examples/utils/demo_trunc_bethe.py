#!/usr/bin/env python
# coding: utf-8

# ## Demo: Memory truncated DMFT simulation

# ### Overview:
# 
# This test program provides a minimal example for the implementation of the truncated Kadanoff Baym Equations.
# 
# We implement a simple Hubbard model on a Bethe lattice and perform an interaction quench from $U=0$ to a finite $U$, using second-order perturbation theory. For small $U$ there is a fast pre-thermalization followed by a slow thermalization dynamics. 
# 
# To be precise, we solve the following self-consistent equations, corresponding to the closed form DMFT selfconsistency equation:
# 
# $$
# G^{-1} = i\partial_t - K,
# $$
# 
# with a self-consistent memory kernel $K(t,\bar{t}) = \Delta(t,\bar{t})+\Sigma(t,\bar{t})$ given by
# 
# $$
# \Delta(t,t^\prime)=G(t,t^\prime), \,\,\,\Sigma(t,t^\prime)=U(t)G(t,t^\prime)G(t,t^\prime)G(t^\prime,t)U(t^\prime)
# $$
# As described for the standard `NESSi`, these equations are solved by a standard forward iteration at each timeslice.  From the converged solution, we can calculate the momentum-dependent Green's function $G_k$, which (due to the $k$-independent self-energy) depends only on the single particle energy $\epsilon_k$, by solving the Dyson equation:
# $$
# G_k^{-1} = i\partial_t - \epsilon_k - \Sigma
# $$
# For the given (particle-hole symmetric) model $\epsilon_k=-2,2$ corresponds to the edge of the nonintreracting bandwidth, while $\epsilon_k=0$ corresponds to the Fermi energy.
# 

# ### Outline of the simulation:
# 
# * To solve the equations with a memory cutoff of `tc` timesteps in the kernel $K$, we will first do a solution of the equations on the full contour up to a given number of timesteps `nt`, where `nt>=tc`. 
# * The result (full Green's functions $G$, $\Sigma$, $G_k$ stored as `cntr::herm_matrix<double>`) is written to an hdf5 file. 
# * To do then a memory truncated evolution from time `tstp=tc+1,...,tmax`, we read in the initial full Green's functions ($G$, $\Sigma$, $G_k$)  from the file, initialize corresponding moving Green's functions windows of type `cntr::herm_matrix_moving<double>`, and do the truncated time evolution

# ### Files
# 
# - `programs/trunc_bethe_start.cpp`: Source code for initial untruncated calculation. GFs stored to hdf5 output  in the working directory.
# - `programs/trunc_bethe.cpp`: Source code for truncated calculation. Reads Green's functions from previous output, and perform the truncated calculation.
# - `utils/demo_trunc_bethe.ipynb`: This notebook.
# - `utils/demo_trunc_bethe.py`: Python script to run the program. Does the same as this Notebook.
# 
# 
# ### Running/Compilation hints
# 
# For the installation, follow the instuction for the installation of the `NESSi` examples. The programs must be compiled with `hdf5=ON`. To run, copy the notebook to a working directory, and adapt the path of the exectutables (see below). The Python script depends on the Python utitilies for reading hdf5 Green's function files. Make sure that `libcntr/python3` is part of the `PYTHONPATH`, , e.g., using the `sys.path.append` code below

# In[ ]:


# -- all Python  modules
import sys
import os
import glob
import numpy as np
import subprocess
# ---- for plotting ----
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap


# In[ ]:


# === Global paths ===
# adapt this as nessesary:
path_to_libcntr = "/path/to/nessi/libcntr/"  
path_to_examples = "/path/to/nessi/examples/"  
# directory to store input and output [relative to directory from which you run this notebook]
output_folder = "./"
os.makedirs(output_folder, exist_ok=True)

# exectutables --- 
nontruncated_simulation = path_to_examples + "exe/trunc_bethe_start.x"
truncated_simulation = path_to_examples + "exe/trunc_bethe.x"

#Import python modules form nessi/
sys.path.append(path_to_libcntr + "python3/") #update PYTHONPATH
from ReadCNTR import *
from ReadCNTRhdf5 import *


# ### Non-truncated simulation: 
# 
# The self-consistent KBE are solved using the standard `NESSI` `dyson` routines, with a accuracy `SolveOrder=MAX_SOLVE_ORDER=5`
# 
# **Parameters:**
# 
# * `U1`: Interaction at $t>0$
# * `beta`: Initial inverse temperature
# * `nt`: Number of timesteps on the inital full KBE simulation.
# * `h`: Timestep
# * `ntau`: Number timesteps on the imaginary axis
# * `BootstrapMaxIter`: Max number of iterations for the bootstrapping phase (time 1,...,SolveOrder)
# * `BootstrapMaxErr`: Error cutoff for the bootstrapping phase
# * `CorrectorSteps`: Fixed Number of selfconsistency iterations at each timstep
# 
# Furthermore, the values $\epsilon_k$ are provided as an `np.array`. The output is stored in a file **[run_name]data.h5**
# 
# **Running the non-truncted simulation**
# 
# The function call is 
# 
# > .../test_truncation_bethe_start.x input-file output-file epsk_1 epsk_2 ...
# 
# where `epsk_1 epsk_2 ...` is a list of the $\epsilon_k$ values for which $G_k$ is computed. If successful, the program generates an hdf5 output file `output-file` which stores the Green's functions. The std-output (convergence error at each timestep and iteration) is redirected to a file `log_file.txt` in the working directory. The function below writes the input and calls the program in this way:

# In[ ]:


def run_start(run_name,nt,ntau,beta,h,U1,CorrectorSteps,BootstrapMaxIter,BootstrapMaxErr,epsk_values):
    # ### Generate the input file in the working directory
    input_file = output_folder + run_name + "input.txt"
    out_file = output_folder + run_name + "data.h5"
    log_file = output_folder + run_name + ".log"
    with open(input_file, 'w') as f:
        f.write(f"__nt= {nt}\n")
        f.write(f"__beta= {beta}\n")
        f.write(f"__h= {h}\n")
        f.write(f"__ntau= {ntau}\n")
        f.write(f"__U1= {U1}\n")
        f.write(f"__CorrectorSteps= {CorrectorSteps}\n")        
        f.write(f"__BootstrapMaxIter= {BootstrapMaxIter}\n")        
        f.write(f"__BootstrapMaxErr= {BootstrapMaxErr}\n")  
    with open(log_file, "w") as log_file_handle:
        subprocess.run([nontruncated_simulation, input_file, out_file, *map(str, epsk_values)],stdout=log_file_handle)
    with open(log_file, "r") as log_file_handle:
        lines = log_file_handle.readlines()
        for line in lines[-10:]:
            print(line, end="")  # avoid double newlines        


# In[ ]:


## run for some params:
run_name="run001"
U1= 1
h= 0.04
BootstrapMaxIter= 100
BootstrapMaxErr= 1e-6
CorrectorSteps= 4
# The values of eps_k for which the simulation is done:
epsk_values = np.array([0.05,0.1,0.15,0.2,0.5,1.0,2.0])  
nk=len(epsk_values)

## params used for the paper (runtime ~ a minute)
nt=800
beta= 100
ntau= 2000

run_start(run_name,nt,ntau,beta,h,U1,CorrectorSteps,BootstrapMaxIter,BootstrapMaxErr,epsk_values)
out_file_start=output_folder + run_name + "data.h5"


# ### Analysis of the decay of the self-energy
# 
# We can now analyze the decay of the self-energy. The following code generates plots of $|\Sigma^{X}(t,t')|$ and $|\Sigma^{X}(t,t')|$ for $X=R,<$  on the last timeslice `nt` (line plots), and in the full $t,t'$ plane (color plots). We extract the self-energy and Green's functions from the hdf5 file using the python helper routines in the modules `ReadCNTR`  and `ReadCNTRhdf5` provided by `NESSi`  (see `NESSi` documentation).
# 
# **Notes:**
# 
# * All functions decay exponentially from the time diagonal, which justifies the truncation approach. 
# * Because the memory integrals are in fact convolutions of the Green's function and the Kernel and not products, one cannot directly estimate the precise necessary memory depth `tc` by looking at the decay of the Kernel. `tc` could be longer and will later be used as a numerical convergence parameter instead.
# * The Kernel from the self-consistent determination of $G$ is $K=\Sigma+G$. Because $\Sigma$ is a pointwise product of Green's functions, $G$ decays more slowly and therefore dominates the decay of the Kernel $K$. 

# In[ ]:


def read_GF_from_datafile(out_file,key):
    with h5py.File(out_file, 'r') as fd:
        # extract the Greenfunction as a list of arrays
        G=read_gf_group(fd[key]) 
        # now G.ret[t1,t2,i,j] contains G^ret(t11,t2)_{ij} for t1>=t2
        # and G.les[t1,t2,i,j] contains G^ret(t1,t2)_{ij} for t2>=t1 (see NESSi documentation)
        # Since here we have scalar 1x1 Greens functions, we extract the 00 component:
        ret=G.ret[:,:,0,0]
        les=G.les[:,:,0,0]
        # symmetrize leser component to fill t1>t2
        lower_triangular_indices = np.tril_indices_from(les, k=-1)
        # Assign the complex conjugate of the upper triangular values to the lower part
        les[lower_triangular_indices] = les.T[lower_triangular_indices].conj()
    return ret,les


# In[ ]:


# --- Custom colormap
colors = [(1, 1, 1), (0.5, 0, 0.5), (1, 0, 0)]  # white → purple → red
custom_cmap = LinearSegmentedColormap.from_list('white_to_bright', colors, N=256)

Sret, Sles = read_GF_from_datafile(out_file_start, 'Sigma')
Gret, Gles = read_GF_from_datafile(out_file_start, 'G')

# ======================  RETARDED COMPONENTS =========================
fig, axes = plt.subplots(1, 2, figsize=(0.8*9, 4), sharex=False, sharey=False)
fig.subplots_adjust(wspace=0.2)  # increase horizontal spacing

# ---- Left panel: |Sigma^R|
abs_Sret = np.abs(Sret)
nt1_L = abs_Sret.shape[0]
im0 = axes[0].imshow(abs_Sret,extent=[0, (nt1_L - 1) * h, 0, (nt1_L - 1) * h],interpolation='none',cmap=custom_cmap,origin='lower')
axes[0].text(0.1, 0.96, r'$|\Sigma^R(t_1,t_2)|$', transform=axes[0].transAxes,fontsize=12, fontweight='bold', va='top')
axes[0].set_xlabel(r"$t_2$")
axes[0].set_ylabel(r"$t_1$")
axes[0].text(0.02, 0.95, "a)", transform=axes[0].transAxes,fontsize=12, fontweight='bold', va='top')

# ---- Right panel: |G^R|
abs_Gret = np.abs(Gret)
nt1_R = abs_Gret.shape[0]
im1 = axes[1].imshow(abs_Gret,extent=[0, (nt1_R - 1) * h, 0, (nt1_R - 1) * h],interpolation='none',cmap=custom_cmap,origin='lower')
axes[1].text(0.1, 0.96, r'$|G^R(t_1,t_2)|$', transform=axes[1].transAxes,fontsize=12, fontweight='bold', va='top')
axes[1].set_xlabel(r"$t_2$")
axes[1].set_ylabel(r"$t_1$")
axes[1].text(0.02, 0.95, "b)", transform=axes[1].transAxes,fontsize=12, fontweight='bold', va='top')

# ---- Explicit colorbars (manual size/position)
fig.canvas.draw()  # update axes positions before reading them

# Left colorbar
cax0 = fig.add_axes([0.4,0.25,0.017,0.2])
fig.colorbar(im0, cax=cax0, orientation='vertical')

# Right colorbar
cax1 = fig.add_axes([0.8,0.25,0.017,0.2])
fig.colorbar(im1, cax=cax1, orientation='vertical')

# ---- Save
png_file = os.path.join(output_folder, "trunc_GR.png")
plt.savefig(png_file, dpi=300, bbox_inches="tight")
plt.show()


# In[ ]:


# --- Custom colormap
colors = [(1, 1, 1), (0.5, 0, 0.5), (1, 0, 0)]  # white → purple → red
custom_cmap = LinearSegmentedColormap.from_list('white_to_bright', colors, N=256)

Sret, Sles = read_GF_from_datafile(out_file_start, 'Sigma')
Gret, Gles = read_GF_from_datafile(out_file_start, 'G')

# ======================  LESSER COMPONENTS =========================
fig, axes = plt.subplots(1, 2, figsize=(0.8*9, 4), sharex=False, sharey=False)
fig.subplots_adjust(wspace=0.2)  # increase horizontal spacing

# ---- Left panel: |Sigma^<|
abs_Sles = np.abs(Sles)
nt1_L = abs_Sles.shape[0]
im0 = axes[0].imshow(abs_Sles,extent=[0, (nt1_L - 1) * h, 0, (nt1_L - 1) * h],interpolation='none',cmap=custom_cmap,origin='lower')
axes[0].text(0.1, 0.96, r'$|\Sigma^<(t_1,t_2)|$', transform=axes[0].transAxes,fontsize=12, fontweight='bold', va='top')
axes[0].set_xlabel(r"$t_2$")
axes[0].set_ylabel(r"$t_1$")
axes[0].text(0.02, 0.95, "a)", transform=axes[0].transAxes,fontsize=12, fontweight='bold', va='top')

# ---- Right panel: |G^<|
abs_Gles = np.abs(Gles)
nt1_R = abs_Gles.shape[0]
im1 = axes[1].imshow(abs_Gles,extent=[0, (nt1_R - 1) * h, 0, (nt1_R - 1) * h],interpolation='none',cmap=custom_cmap,origin='lower')
axes[1].text(0.1, 0.96, r'$|G^<(t_1,t_2)|$', transform=axes[1].transAxes,fontsize=12, fontweight='bold', va='top')
axes[1].set_xlabel(r"$t_2$")
axes[1].set_ylabel(r"$t_1$")
axes[1].text(0.02, 0.95, "b)", transform=axes[1].transAxes,fontsize=12, fontweight='bold', va='top')

# ---- Explicit colorbars (manual size/position)
fig.canvas.draw()  # update axes positions before reading them

# Left colorbar
cax0 = fig.add_axes([0.4,0.25,0.017,0.2])
fig.colorbar(im0, cax=cax0, orientation='vertical')

# Right colorbar
cax1 = fig.add_axes([0.8,0.25,0.017,0.2])
fig.colorbar(im1, cax=cax1, orientation='vertical')

# ---- Save
png_file = os.path.join(output_folder, "trunc_Gles.png")
plt.savefig(png_file, dpi=300, bbox_inches="tight")
plt.show()


# ### Analysis of the decay of the Greens functions $G_k$
# 
# 
# It is important to note that for the truncation to work, it is sufficient that the Kernel decays. The Green's function which is a solution of the Dyson equation does not need to decay quickly. This is evident for the momentum resolved Green's functions, for which the kernel $\Sigma$ has a quick decay but the Green's functions can decay only very slowly, in particular at the Fermi energy.  The following code generates plots of $|G_k^{R,<}(t,t')|$, $|G^{R,<}(t,t')|$, and $|\Sigma^{R,<}(t,t')|$  on the last timeslice `nt` (line plots).

# In[ ]:


fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharex=False, sharey=False)

## line plot of G and Sigma, (from previously read matrices)
nt1 = Sret.shape[0]
j = nt1 - 1  # last index
# Build x and y data
i_vals = np.arange(nt1)
x = j-i_vals 

axes[0].set_yscale('log')
axes[0].plot(x*h, np.abs(Sret[j,:]), 'b-', label=r"$|{\Sigma}^R(t, t')|$")
axes[0].plot(x*h, np.abs(Sles[j,:]), 'b--', label=r"$|{\Sigma}^<(t, t')|$")
axes[0].plot(x*h, np.abs(Gret[j,:]), 'r-', label=r"$|G^R(t, t')|$")
axes[0].plot(x*h, np.abs(Gles[j,:]), 'r--',label=r"$|G^<(t, t')|$")
axes[0].set_xlabel(r"$t - t'$")
axes[0].set_ylabel(" ")
axes[0].text(0.2, 0.95, "a)", transform=axes[0].transAxes,fontsize=12, fontweight='bold', va='top')
#axes[0].grid(True, which='both', linestyle='--', alpha=0.5)
axes[0].legend()

# ========== plot of momentum dependent Green's functions

cmap = matplotlib.colormaps['coolwarm']

# Dictionary to store k-color mapping
k_color_map = {}

# loop over k values, extract data from HDF5 file:
for k,eps_k in enumerate(epsk_values):
    label = f"lattice_functions/{k}/Gk"
    Gret_k, _ = read_GF_from_datafile(out_file_start,label)  
    nt1 = Gret_k.shape[0]
    j = nt1 - 1
    t_diff = (j - np.arange(nt1)) * h  # x-axis, same for all
    y = np.abs(Gret_k[j, :])
    axes[1].plot(t_diff, y,color=cmap(k/nk), linewidth=2,label=rf"$\epsilon_k = {eps_k:.2f}$")

axes[1].legend()
axes[1].set_yscale('log')
axes[1].set_xlabel(r"$t - t'$")
axes[1].set_ylabel(r"$|G^R_k(t, t')|$")
axes[1].text(0.1, 0.8, "b)", transform=axes[1].transAxes,fontsize=12, fontweight='bold', va='top')

plt.tight_layout()

# ---- Save
png_file = os.path.join(output_folder, "trunc_Gk.png")
plt.savefig(png_file, dpi=300, bbox_inches="tight")

plt.show()



# ### The truncated simulation
# 
# **Parameters:**
# 
# * `CorrectorSteps`: Fixed Number of selfconsistency iterations at each timestep
# * `tc`: The truncation time. `tc<=nt` required, such that the code can read the inital Green's function from file.
# * `tmax`: Simulation is performed for times `tc+1,...,tmax`
# * `out_every`: frequency at which timeslices are written to file
# * `U1`, `h`, as well as the $k$ values are read from the output  **[run_name]data.h5** of the startup
# 
# 
# The function call is 
# 
# > .../test_truncation_bethe.x input-file data_file nk_out_file
# 
# `data_file` is the hdf5 output generated by the startup run.
# 
# If successful, the program generates an text output file `out_file_nk` which stores the momentum occupations, 
# 
# $$
# \rho_k(t) = -i G_k^<(t,t)
# $$
# 
# in the text form with columns $\epsilon_k$, timestep, $\rho_k(t)$, and an HDF5 file `out_file_slices` which containes selected timeslices `tc, tc+out_every`, `tc, tc+2*out_every`,...

# In[ ]:


def trunctated_simulation(run_name,data_file,out_file_nk,out_file_slices,tc,tmax,out_every,CorrectorSteps):
    input_file = output_folder + run_name + "input.txt"
    #start_file = output_folder + run_name + "data.h5"
    #out_file = output_folder + run_name + nk_out_file
    log_file = output_folder + run_name + "trunc.log"
    #generate iuput 
    with open(input_file, 'w') as f:
        f.write(f"__tmax= {tmax}\n")
        f.write(f"__tc= {tc}\n")
        f.write(f"__out_every= {out_every}\n")
        f.write(f"__CorrectorSteps= {CorrectorSteps}\n")
    with open(log_file, "w") as log_file_handle:
        subprocess.run(
            [truncated_simulation, input_file, data_file, out_file_nk, out_file_slices],
            stdout=log_file_handle
        )


# In[ ]:


# function to parse the nk-output
def read_nk_data(nk_out_file, epsk_values, nt):
    nk = len(epsk_values)
    rhok = np.zeros((nk, nt + 1))
    with open(nk_out_file, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) != 3:
                continue  # skip malformed or comment lines
            k_val, t, nk_val = parts
            k_val = float(k_val)
            t = int(t)
            nk_val = float(nk_val)
            # Use np.isclose in case of float precision mismatch
            k_idx = np.where(np.isclose(epsk_values, k_val))[0][0]
            rhok[k_idx, t] = nk_val
    return rhok, epsk_values



# In[ ]:


# run the simulation for three tc's can take few minutes; for quick check, reduce tmax

tmax=30000 #used for the plots in the manuscript
tc_vals=np.array([200,400,600,800]) 
CorrectorSteps= 4
out_every=tmax+1  # no intermediade output

for tc in tc_vals:
    out_file_nk = output_folder + run_name + f"tc{tc}.txt"
    out_file_slices = output_folder + run_name + f"slices.h5" # will not be used here
    data_file = output_folder + run_name + "data.h5"   
    print(f"simulation with tc={tc}")
    trunctated_simulation(run_name,data_file,out_file_nk,out_file_slices,tc,tmax,out_every,CorrectorSteps)


# ### Plot of the momentum occupation:
# 
# After running the simulation, this block generates a plot of $\rho_k(t)$ for various $k$ values. One can see a two-stage dynamics, with a prethermalization after times of order $1$, and a later thermalization. If initial temperature $T=1/\beta$ is sufficiently low, one can see that in the prethermal states the momentum distribution has a pronounced step at the Fermi energy, by comparing $\rho_k$ at the smallest $\epsilon_k=0.05$ to $\rho_{k_f}=0.5$ at the Fermi energy $\epsilon_{k_F}$. 
# 
# By computimg the data with various `tc`, one can see a convergence with `tc` much smaller than the maximum simulation time `tmax`. (Note the logarithmic time axis).
# 

# In[ ]:


# Define line styles for each tc value
line_styles = [':', '--', '-', ':']

cmap = matplotlib.colormaps['coolwarm']


# Create a figure for plotting
plt.figure(figsize=(9, 4.5))

# Dictionary to store k-color mapping
k_color_map = {}

# Loop over all tc values and plot
for tc_idx, tc in enumerate(tc_vals[:3]):
    nk_out_file = output_folder + run_name + f"tc{tc}.txt"

    # Read data for the current tc value
    nkvals, k_values = read_nk_data(nk_out_file, epsk_values, tmax)

    # Loop over k values and plot the corresponding nk[k][t] for each
    for idx, k in enumerate(k_values):
        # Plot with the same color but different line style for each tc
        line, = plt.plot(np.arange(tmax + 1) * h, nkvals[idx], 
                         color=cmap((1.0*idx)/nk),
                         linestyle=line_styles[tc_idx % len(line_styles)], 
                         label=f'{k} (tc={tc})')

        # Store the color used for this k value (only once)
        if k not in k_color_map:
            k_color_map[k] = line.get_color()

# Set y-axis range
plt.ylim(-0.01, 0.55)
plt.xscale('log')

# Labels and title
plt.xlabel('t')
plt.ylabel(r'$\rho_k(t)$')
#plt.title(r'\rho_k(t)$')

# **Legend for line styles (representing tc values)**
line_types_legend = [Line2D([0], [0], color='black', lw=2, linestyle=style) for style in line_styles[:3]]
line_labels = [fr'$t_c={tc * h}$' for tc in tc_vals[:3]]

# **Legend for k values (representing color code)**
k_legend_handles = [Line2D([0], [0], color=color, lw=2) for color in k_color_map.values()]
k_legend_labels = [fr'$\epsilon_k={k}$' for k in k_color_map.keys()]

# Add legends to the plot
legend1 = plt.legend(handles=line_types_legend, labels=line_labels, loc='upper left', bbox_to_anchor=(0.02, 1.0))
legend2 = plt.legend(handles=k_legend_handles, labels=k_legend_labels, loc='upper left', bbox_to_anchor=(0.2, 1.0))

# Add the first legend back so both are visible
plt.gca().add_artist(legend1)

#plt.grid(True)

# ---- Save
png_file = os.path.join(output_folder, "trunc_convergence.png")
plt.savefig(png_file, dpi=300, bbox_inches="tight")


plt.show()


# ### A quick demo of the HDFT output at slices:
# 
# We run the simulation from above, but write a timeslice every 1000 steps.
# 
# We can then extract the slices using the `read_herm_matrix_timestep_moving_moving_group` 
# function from `readCNTRhdf5` to read the slices. Not much is happening, as local Green's 
# functions vary weakly after a quench. This just demonstrates the i/o.

# In[ ]:


tmax=10000 #used for the plots in the manuscript
tc=200
out_every=1000
CorrectorSteps= 4
out_file_nk = output_folder + run_name + "nk.txt"
out_file_slices = output_folder + run_name + f"slices.h5" 
data_file = output_folder + run_name + "data.h5"   
trunctated_simulation(run_name,data_file,out_file_nk,out_file_slices,tc,tmax,out_every,CorrectorSteps)


# In[ ]:


Gret_slice=[]
Gles_slice=[]
ts=np.arange(tc,tmax,out_every)

with h5py.File(out_file_slices, 'r') as fd:
    for t in ts:
        Gd=read_herm_matrix_timestep_moving_group(fd[f"t{t}/G"])
        # Since here we have scalar 1x1 Greens functions, we extract the 00 component:
        Gret_slice.append(Gd.ret[:,0,0])
        Gles_slice.append(Gd.les[:,0,0])

## ==== plot

fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharex=False, sharey=False)

cmap = matplotlib.colormaps['coolwarm']

tvalues=np.arange(tc+1)*h

for idx,t in enumerate(ts):
    axes[0].plot(tvalues, np.imag(Gles_slice[idx]) , '-', color=cmap((t-tc)/(tmax-tc)), label=rf"$t = {t*h:.2f}$")
    axes[1].plot(tvalues, np.imag(Gles_slice[idx]-Gles_slice[0]) , '-', color=cmap((t-tc)/(tmax-tc)), label=rf"$t = {t*h:.2f}$")

axes[0].set_xlabel(r"$s$")
axes[0].set_ylabel(r"Im$G^<(t,t-s)$")
axes[1].set_xlabel(r"$s$")
axes[1].set_ylabel(r"Im$G^<(t,t-s)-G^<(t_c,t_c-s)$")

axes[0].legend()
plt.tight_layout()
plt.show()




# In[ ]:


# clear all data files produced by this worksheet:
os.system(f"rm {output_folder}run001* ")


# In[ ]:




