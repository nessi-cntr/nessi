libcntr library
===============

This is the **libcntr** library, the core part of the **NESSi** software package. The full documentation and detailed build instructions
can be found at the [official webpage](http://www.nessi.tuxfamily.org).

Build instructions
==================

### Prerequisites
* [cmake](https://cmake.org)
* [eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) 
* [hdf5](https://www.hdfgroup.org/solutions/hdf5/) - optional, but recommended

### Configuration and installation

We use the [cmake](https://cmake.org) build system to generate Makefiles. We recommend creating a configuration script 
with all options:

```
# From the root directory .../libcntr/
# Create configure script
vim configure.sh
```

The basic configure script has the following structure:

```
CC=[C compiler] CXX=[C++ compiler]
cmake \
 -DCMAKE_INSTALL_PREFIX=[install directory] \
 -DCMAKE_BUILD_TYPE=[Debug|Release] \
 -DCMAKE_INCLUDE_PATH=[include directory] \
 -DCMAKE_LIBRARY_PATH=[library directory] \
 -DCMAKE_CXX_FLAGS="[compiling flags]" \
 ..
```

In the first line, the C and C++ compiler are set. We have tested the full library for both the GNU compilers (`CC=gcc CXX=g++`) as well as the Intel compilers (`CC=icc CXX=icpc`). The install directory (for instance `/home/opt`) is defined by the CMake variable `CMAKE_INSTALL_PREFIX`. Debugging tools are switched on by setting `CMAKE_BUILD_TYPE=Debug`; otherwise, all assertions and sanity checks are turned off. The code is significantly faster in release mode (`CMAKE_BUILD_TYPE=Release`) and thus recommended for a production run. The debug mode, on the other hand, turns on assertions (implemented as C++ standard assertions) checking a consistency of the input for all major routines. 
The path to the libraries that **libcntr** depends upon ([eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page)  and, optionally, [hdf5](https://www.hdfgroup.org/solutions/hdf5/)) are provided by specifying the include directory `CMAKE_INCLUDE_PATH` and the library path `CMAKE_LIBRARY_PATH`. Finally, the compilation flags are specified by `CMAKE_CXX_FLAGS`. To compile libcntr, the flags should include

```
-std=c++11 -fpermissive
```

if the GNU C++ compiler is used. No special flags are required when using the Intel C++ compiler.

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

After successful configuration (which generates the make files), compile the library by

```
make
```


HDF5
====

To enable hdf5 support add the `hdf5` option to the cmake configure script:

    -Dhdf5=ON

MPI and OpenMP
==============

To turn on OpenMP and/or MPI parallelization define the options `omp` and/or `mpi` in the cmake step, respectively, and specify your C and C++ MPI compilers by `CC=mpicc CXX=mpix++` (or similar, depending on your system). Add the `omp` and `mpi` option to the configure script:

    -Domp=ON \
    -Dmpi=ON


doxygen
==============

The **libcntr** is fully documented using the automated documentation tool [doxygen](http://www.doxygen.nl). To build the documentation,
add the cmake variable

    -DBUILD_DOC=ON

to the configure script. Upon running `make`, the documentation will be generated under `doc/html/index.html`.


Tests
=====

We also provide a test suite for checking the functionality of every major routine in **libcntr** based on the Catch library. For running the tests, simply run


    make test

  
After completing all test, the message `All tests` passed indicates that the compiled version of **libcntr** is fully functional. If compiled with MPI support, the MPI-based functions can be tested by running


    make test_mpi

Example scripts
===============

On MacOSX and using the GNU compilers, an example configure script with all options would look like the following:

```
CC=mpicc CXX=mpicxx \
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/opt \
    -DCMAKE_INSTALL_NAME_DIR=$HOME/opt/lib \
    -DCMAKE_MACOSX_RPATH=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_DOC=ON \
    -Domp=ON \
    -Dhdf5=ON \
    -Dmpi=ON \
    -DCMAKE_INCLUDE_PATH=/opt/local/include \
    -DCMAKE_LIBRARY_PATH=/opt/local/lib \
    -DCMAKE_CXX_FLAGS="-std=c++11 -fpermissive" \
    ..
```

The dynamic linking scheme under MacOSX requires specifying `-DCMAKE_INSTALL_NAME_DIR=$HOME/opt/lib` as the directory where the library `libcntr.dylib` will be installed. In this example, all libraries (eigen3, hdf5) have been installed under the prefix `opt/local`.

Under Linux, the above job script would look like

```
CC=mpicc CXX=mpicxx \
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/opt \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_DOC=ON \
    -Domp=ON \
    -Dhdf5=ON \
    -Dmpi=ON \
    -DCMAKE_INCLUDE_PATH=/usr/local/include \
    -DCMAKE_LIBRARY_PATH=/usr/local/lib \
    -DCMAKE_CXX_FLAGS="-std=c++11 -fpermissive" \
    ..
```

Here we have assumed that the dependencies are installed under `/usr/local/`.


