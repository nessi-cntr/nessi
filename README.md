# NESSi v2.0 - The NonEquilibrium Systems Simulation package!
## What is NESSi? 
`NESSi` is an open-source software package for the manipulation of nonequilibrium Green’s functions defined on the Kadanoff-Baym contour. The Green's function method in its time-dependent formulation is a versatile framework for the solution of interacting many-body problems out of equilibrium. `NESSi` provides classes representing the various types of Green’s functions, implements the basic operations on these functions and allows to solve the corresponding equations of motion. The library is aimed at the study of transient dynamics from an initial equilibrium state, induced by time-dependent model parameters.    
### Overview
- `NESSi` provides tools for constructing Feynman diagram and solving equations of motion for non-equilibrium Green's functions on the Kadanoff-Baym contour
- `NESSi` is based on high-order quadrature rules: for $N$ time slices, the error scales like up to $\mathcal{O}(N^{-p})$ with $p$ up to $7$.
- Efficient distributed-memory parallelization over reciprocal space allows large-scale calculations on extended systems.
- The first extension of `NESSi` enables memory truncated time propagation and steady-state calculations

## Core functionalities
The `libcntr` library provides highly accurate methods for calculating nonequilibrium Green's functions and more. A brief overview of the core routines is presented below. All routines work for fermions and bosons.
1. Construct the free Green's functions for a general time-dependent Hamiltonian.
2. Solve the Dyson equation along the full Kadanoff-Baym contour $`{\cal C}`$: $` [i\partial_t - h(t)] G (t,t') - [\Sigma * G](t,t') = \delta_{\cal C} (t,t') `$ for a given self-energy $\Sigma(t,t')$, also in a memory truncated propagation scheme. In particular, thermal equilibrium and time evolution are treated on equal footing.
3. Solve the contour integral equations of the type $` G(t,t^\prime) + [F\ast G](t,t^{\prime}) = Q(t,t^\prime)`$ for given kernel $`F(t,t^\prime)`$ and $`Q(t,t^\prime)`$. A typical example is the self-consistent $`GW`$ approximation, where the screened interaction obeys the Dyson equation $`W = V + V\ast \Pi \ast W`$. Here $V$ denotes the bare Coulomb interaction, while $\Pi$ stands for the irreducible polarization.
4. Compute the convolution $`[A\ast B](t,t^\prime) = \int_{\mathcal{C}} d\bar{t} \, A(t,\bar{t}) B(\bar{t},t^\prime)`$ of two contour functions $`A(t,t^\prime)`$ and $`B(t,t^\prime)`$, which is an essential part of solving of the integro-differential equations.
5. Construct the bubble diagrams of the type $` C(t,t^\prime)= i A(t,t^\prime) B(t^\prime, t)`$ or $`C(t,t^\prime)= i A(t,t^\prime) B(t, t^\prime)`$.

`NESSi` version 2.0 contains the memory truncated time propagation and the steady state formalism. This update extends `NESSi` by the following core functionalities:
1. Solves the truncated Dyson equation along the full Kadanoff-Baym contour $`{\cal C}`$: $` [i\partial_t - h(t)] G (t,t') - [\Sigma * G](t,t') = \delta_{\cal C} (t,t') `$ for a given self-energy $`\Sigma(t,t')`$.
2.  Solves the truncated contour integral equations of the type $` G(t,t^\prime) + [F\ast G](t,t^\prime) = Q(t,t^\prime) `$ for given kernel $`F(t,t^\prime)`$ and $`Q(t,t^\prime)`$.
3. Construct the bubble diagrams of the type $` C(t,t^\prime)= i A(t,t^\prime) B(t^\prime, t) \quad \text{or}\quad C(t,t^\prime)= i A(t,t^\prime) B(t, t^\prime)`$ for truncated Green's functions.
4. Solve the Dyson equation along the full Kadanoff-Baym contour $`{\cal C}`$: $` [i\partial_t - h(t)] G (t,t') - [\Sigma * G](t,t') = \delta_{\cal C} (t,t') `$ in the steady state for a given self-energy $`\Sigma(t,t')`$. The routine evaluates the Dyson equation in frequency domain and calcualtes the necessary Fourier integral transforms.
5. Construct the equilibrium steady-state Green's function $` G(t) = \int d \omega \,A(\omega) g_{\omega}(t) `$ for a given Density of states $` A(\omega) `$.
6. Construct the bubble diagrams of the type $` C(t,t^\prime)= i A(t,t^\prime) B(t^\prime, t) \quad \text{or}\quad C(t,t^\prime)= i A(t,t^\prime) B(t, t^\prime)`$ for steady-state Green's functions.
## How to cite
Please cite the following paper whenever you use parts of NESSi:

*M. Schüler, D. Golež, Y. Murakami, N. Bittner, A. Herrmann, H. U. R. Strand, P. Werner, and M. Eckstein, NESSi: The Non-Equilibrium Systems Simulation package, Comput. Phys. Commun. 257, 107484 (2020)*

Please cite in addition the following paper whenever you use truncated or steady-state objects from the `NESSi` package:

Version 2.0.0 paper, to be published

## Documentation
- Documentation and manual available on [NESSi website](https://nessi.readthedocs.io/en/latest/)
- Documentation and manual can be [locally build](#Manual-and-documentation)


## Installation
To install the `NESSi` library, follow the [instructions](#NESSi) below.


### NESSi
The `nessi` directory contains two folders: `libcntr`, which contains the library along with a collection of python tools, and `examples`, which includes example programs.

The library `libcntr` depends on the following three libraries: [eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) (required), [HDF5](https://www.hdfgroup.org/solutions/hdf5/)(optional, recommended),  [FFTW3](https://www.fftw.org)(optional, required for steady state code).

For compiling and installing the `libcntr` library, we use the `CMake` building environment to generate system specific make files. Several `CMake` variables need to be set to install `libcntr`; their are passed as `-Dvar=[value]`. `CMake` can be called directly from the terminal; however, it is more convenient to create a configure script with all variables and compile options. We suggest a script of the following structure:

    CC=[C compiler] CXX=[C++ compiler]
    cmake \
     -DCMAKE_INSTALL_PREFIX=[install directory] \
     -DCMAKE_BUILD_TYPE=[Debug|Release] \
     -Domp=[ON|OFF] \
     -Dhdf5=[ON|OFF] \
     -Dmpi=[ON|OFF] \
     -DBUILD_DOC=[ON|OFF] \
     -Dness=[ON|OFF] \
     -DCMAKE_INCLUDE_PATH=[include directory] \
     -DCMAKE_LIBRARY_PATH=[library directory] \
     -DCMAKE_CXX_FLAGS="[compiling flags]" \
     ..

, which in the following will be refered to as `configure.sh`.
In the first line, the C and C++ compiler are set. We have tested the full library for both the GNU compilers (` CC=gcc CXX=g++ `) as well as the Intel compilers (`CC=icc CXX=icpc`). The install directory (for instance `/home/opt`) is defined by the `CMake` variable `CMAKE_INSTALL_PREFIX`. Debugging tools are switched on by setting `CMAKE_BUILD_TYPE=Debug`; otherwise, all assertions and sanity checks are turned off. The code is significantly faster in release mode (`CMAKE_BUILD_TYPE=Release`) and thus recommended for a production run. The debug mode, on the other hand, turns on assertions (implemented as C++ standard assertions) checking a consistency of the input for all major routines. The following three lines trigger optional (but recommended) functionalities. Setting `omp=ON` turns on the compilation of routines parallelised with `openMP`, while setting `mpi=ON` is required for compiling distributed-memory routines based on MPI. In this case, MPI compilers have to be specified in the first line.  `Dness=ON` turns on the compilation of the steady-state routines, if the library is build with `omp=ON` and `ness=ON`, parallel routines in the steady state code are installed. In addition it requires the FFTW library to be compiled with `--enable-threads` (see [FFTW](https://www.fftw.org)). Finally, `hdf5=ON` activates the usage of the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) library. Here, we recommend to install a version, which is no older than 1.12. The path to the libraries that `libcntr` depends upon ([eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) and, optionally, [HDF5](https://www.hdfgroup.org/solutions/hdf5/)) are provided by specifying the include directory `CMAKE_INCLUDE_PATH` and the library path `CMAKE_LIBRARY_PATH` (note that for using `openMP` the respective paths also need to be included). Finally, the compilation flags are specified by `CMAKE_CXX_FLAGS`. To compile `libcntr`, the flags should include

    -std=c++11
    
As the next step, create a build directory (for instance cbuild). Navigate to this directory and run the configure script:

    sh ../configure.sh
    
After successful configuration (which generates the make files), compile the library by

    make
    
and install it to the install directory by

    make install
    
We have tested the installation procedure on the following platforms: MacOSX with [MacPorts](https://www.macports.org), MacOSX with [homebrew](https://brew.sh), and Linux (Ubuntu, Debian, CentOS and Arch Linux) using the respective native package manager. Detailed installation instructions (with examples) can be found below.

#### On MacOSX with MacPorts

We recommend using the `gcc` compiler. As an example we show the installation with `gcc9` and `openmpi`. The dependencies of `libcntr` are installed via

    sudo port install gcc9 eigen3-devel hdf5
    sudo port install openmpi-devel-gcc9 

For installing doxygen to be able to build the documentation, run

    sudo port install doxygen graphviz

The `python` packages `numpy`, `scipy`, `matplotlib` and `h5py` can be installed using the system's `python` package manager `pip` or using `MacPorts`.

Below we provide an example configure script based on default path (`/opt/local/`) of `MacPorts`:

    CC=mpicc CXX=mpicxx
    cmake \
     -DCMAKE_INSTALL_PREFIX=$HOME/opt \
     -DCMAKE_INSTALL_NAME_DIR=$HOME/opt/lib \
     -DCMAKE_BUILD_TYPE=Release \
     -Domp=ON \
     -Dhdf5=ON \
     -Dmpi=ON \
     -DBUILD_DOC=ON \
     -Dness=ON \
     -DCMAKE_INCLUDE_PATH=/opt/local/include \
     -DCMAKE_LIBRARY_PATH=/opt/local/lib \
     -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
     ..

Note that under MacOSX the installation directory (`$HOME/opt`) has to be provided explitly for other programs to find the shared library `libcntr`, which is accomplished by setting `CMAKE_INSTALL_NAME_DIR=$HOME/opt/lib`.

We recommend using `MacPorts` as it currently supports all compiling options directly.

#### On MacOSX with homebrew

We have tested the compilation based Apple `clang`. Due to some legacy dependence, run

    brew unlink python@2

first. The dependencies are installed by

    brew install eigen hdf5 open-mpi

and

    brew install doxygen graphviz

for building the documentation.

Using `gcc` is also possible, but `openmpi` has to be build from source in this case.

The following example configure script illustrates how to install `libcntr`:

    CC=mpicc CXX=mpicxx
    cmake \
     -DCMAKE_INSTALL_PREFIX=$HOME/opt \
     -DCMAKE_INSTALL_NAME_DIR=$HOME/opt/lib \
     -DCMAKE_BUILD_TYPE=Release \
     -Domp=OFF \
     -Dhdf5=ON \
     -Dmpi=ON \
     -DBUILD_DOC=ON \
     -Dness=ON \
     -DCMAKE_INCLUDE_PATH=/usr/local/include \
     -DCMAKE_LIBRARY_PATH=/usr/local/lib \
     -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
     ..

Note that we disabled `openMP` support by `omp=OFF`. This is due to native Apple `clang` compiler not directly supporting `openMP`. Specific work-arounds exist, see for instance [here](https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp).

#### On Linux

The installation works similar under linux. For Debian or Ubuntu the dependencies can be installed by

    apt-get install -y --allow-unauthenticated  libhdf5-serial-dev libopenmpi-dev libeigen3-dev doxygen graphviz cmake

while

    pacman -Sy --noconfirm  hdf5 gcc openmpi eigen doxygen graphviz make cmake

on Arch Linux. For RedHat-based systems like CentOS use `yum` to install the required packages. The include and library path are usually automatically set, such that the following example configure script would work

    CC=mpicc CXX=mpicxx
    cmake \
     -DCMAKE_INSTALL_PREFIX=$HOME/opt \
     -DCMAKE_BUILD_TYPE=Release \
     -Domp=ON \
     -Dhdf5=ON \
     -Dmpi=ON \
     -DBUILD_DOC=ON \
     -Dness=ON \
     -DCMAKE_INCLUDE_PATH="" \
     -DCMAKE_LIBRARY_PATH="" \
     -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
     ..

Note that setting `CMAKE_INSTALL_NAME_DIR` is not required if the installation path (e. g. `$HOME/opt/lib`) can be found under `LD_LIBRARY_PATH`.

We have further compiled and installed the code on compute clusters with CentOS, where a module system provided the required libraries. We have also tested the compilation with Intel compilers by setting `CC=mpiicc CXX=mpiicpc`.

#### Creating a custom program
As a minimal working example, we provide a program which creates a Greens function and prints out its properties. First, we define Greens function `G` with given properties. After the Greens function object is created the program prints out corresponding attributes. Only one header (`cntr.hpp`) needs to be included to make all subroutines in namespace `cntr::` available.
    
    #include <sys/stat.h>
    #include <iostream>
    #include "cntr/cntr.hpp"

    int main(int argc,char *argv[]){
        GREEN G;
        int Nt=100;   // # of timesteps
        int Ntau=100; // number of imaginary time steps
        int size=3;   // Size of Greens function

        G=GREEN(Nt,Ntau,size,FERMION);

        std::cout << "Number of timesteps: " << G.nt() << std::endl;
        std::cout << "Number of Matsubara points: " << G.ntau() << std::endl;
        std::cout << "Size of Greens function: " << G.size1() << std::endl;
    return 0;
    }

#### Test suite and documentation

We also provide a test suite for checking the functionality of every major routine in `libcntr` based on the Catch library. For running the tests, simply run

    make test

After completing all test, the message

    All tests passed
    
indicates that the compiled version of `libcntr` is fully functional. If compiled with MPI support, the MPI-based functions can be tested by running

    make test_mpi

The C++ code is fully documented using the automatic documentation tool [doxygen](http://www.doxygen.nl). For generating the documentation, specify `BUILD_DOC=ON` in the configure script. The documentation is then generated automatically along the compilation. After a successful build, the html documentation can be found in `nessi/libcntr/doc/html/index.html`.


### Manual and documentation
To generate the manual locally, follow the instructions in the `nessi-cntr/nessi_doc` repository, and additionally the local [doxygen](http://www.doxygen.nl) documentation of the code (after compiling with `BUILD_DOC=ON`) can be found in `nessi/libcntr/doc/html/index.html`.


