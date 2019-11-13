# Example programs

This directory contains a number of example programs, demonstrating the usage of the **libcntr** library and tools.
A detailed discussion is presented in [Michael Schüler, Denis Golež, Yuta Murakami, Nikolaj Bittner, Andreas Herrmann, Hugo U. R. Strand, Philipp Werner, Martin Eckstein, arXiv:1911.01211 [cs]](http://arxiv.org/abs/1911.01211)
A short introduction and instructions to run the programs can be found on the [official webpage](http://www.nessi.tuxfamily.org).

Content
==================



Build instructions
==================

### Prerequisites
* [cmake](https://cmake.org)
* [eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) 
* [hdf5](https://www.hdfgroup.org/solutions/hdf5/) - optional, but recommended
* MPI C++ compiler - for MPI example
* [libcntr](libcntr/)

The following table summarizes which programs depend on the hdf5 library or MPI:

Program | Dependencies
------------ | -------------
hubbard_chain_2b.x | -  
hubbard_chain_gw.x | - 
hubbard_chain_tpp.x | -
test_equilibrium.x | -
test_nonequilibrium.x | -
integration.x | -
Holstein_bethe_Nambu_Migdal.x | hdf5
Holstein_bethe_Nambu_uMig.x | hdf5
Holstein_bethe_Migdal.x | hdf5
Holstein_bethe_uMig.x | hdf5
Holstein_impurity_singlebath_Migdal.x | hdf5
Holstein_impurity_singlebath_uMig.x | hdf5
gw.x | hdf5, MPI


### Configuration and installation

As for the [libcntr](libcntr/) library, we use the [cmake](https://cmake.org) build system. 
We recommend a configure script (`examples/configure.sh`) similar to the following:

```
CC=[C compiler] CXX=[C++ compiler]
cmake \
 -DCMAKE_BUILD_TYPE=[Debug|Release] \
 -Domp=[ON|OFF] \
 -Dhdf5=[ON|OFF] \ 
 -Dmpi=[ON|OFF] \
 -DCMAKE_INCLUDE_PATH=[include directory] \
 -DCMAKE_LIBRARY_PATH=[library directory] \
 -DCMAKE_CXX_FLAGS="[compiling flags]" \
 ..
```

The options and compiler flags are similar to the compilation of [libcntr](libcntr/) . In addition to the include and 
library paths provided for the installation of **libcntr**, the path to the library has to be provided. For example, if the eigen3 and the hdf5 library are installed under `/usr/local/`, while **libcntr** was installed under `$HOME/opt/`, one would specify
```
-DCMAKE_INCLUDE_PATH="/usr/local/include;$HOME/opt/include" \
-DCMAKE_LIBRARY_PATH="/usr/local/lib;$HOME/opt/lib" \
```

As the next step create a build directory

```
# From the root directory .../libcntr/
# Create build directory
mkdir cbuild
cd cbuild
```

and run the configure script:

```
sh ../configure.sh
```

After successful configuration, the programs are compiled by

```
make
```

The executables are located in `examples/exe`.